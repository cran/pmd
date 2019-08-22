## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----demodata------------------------------------------------------------
library(pmd)
data("spmeinvivo")
str(spmeinvivo)

## ----rtg-----------------------------------------------------------------
pmd <- getpaired(spmeinvivo, rtcutoff = 10, ng = 10)
plotrtg(pmd)

## ----pmd-----------------------------------------------------------------
plotpaired(pmd)

## ----pmdindex------------------------------------------------------------
# show the unique PMD found by getpaired function
for(i in 1:length(unique(pmd$paired$diff2))){
        diff <- unique(pmd$paired$diff2)[i]
        index <- pmd$paired$diff2 == diff
        plotpaired(pmd,index)
}

## ----std-----------------------------------------------------------------
std <- getstd(pmd)

## ----stdplot-------------------------------------------------------------
plotstd(std)

## ----stdrtplot-----------------------------------------------------------
par(mfrow = c(2,3))
plotstdrt(std,rtcluster = 23,main = 'Retention time group 23')
plotstdrt(std,rtcluster = 9,main = 'Retention time group 9')
plotstdrt(std,rtcluster = 18,main = 'Retention time group 18')
plotstdrt(std,rtcluster = 67,main = 'Retention time group 67')
plotstdrt(std,rtcluster = 49,main = 'Retention time group 49')
plotstdrt(std,rtcluster = 6,main = 'Retention time group 6')

## ----gettarget-----------------------------------------------------------
# you need retention time for independent peaks
index <- gettarget(std$rt[std$stdmassindex])
# output the ions for each injection
table(index)
# show the ions for the first injection
std$mz[index==1]
std$rt[index==1]

## ----pca-----------------------------------------------------------------
library(enviGCMS)
par(mfrow = c(1,2),mar = c(4,4,2,1)+0.1)
plotpca(std$data,lv = as.numeric(as.factor(std$group)),main = substitute(paste(italic('in vivo'), " SPME samples(all peaks)")))
plotpca(std$data[std$stdmassindex,],lv = as.numeric(as.factor(std$group)),main = substitute(paste(italic('in vivo'), " SPME samples(selected peaks)")))

## ----comp----------------------------------------------------------------
stdcluster <- getcluster(std)
# extract pseudospectra for std peak 71
plot(stdcluster$cluster$mz[stdcluster$cluster$i==71],stdcluster$cluster$ins[stdcluster$cluster$i==71],type = 'h',xlab = 'm/z',ylab = 'intensity',main = 'pseudospectra for GlobalStd peak 711')
# export peaks with the highest intensities in each GlobalStd peaks groups.
data <- stdcluster$data[stdcluster$stdmassindex2,]

## ----corcomp-------------------------------------------------------------
corcluster <- getcorcluster(spmeinvivo)
par(mfrow = c(1,3),mar = c(4,4,2,1)+0.1)
plotpca(std$data,lv = as.numeric(as.factor(std$group)),main = substitute(paste(italic('in vivo'), " SPME samples(all peaks)")))
plotpca(std$data[std$stdmassindex,],lv = as.numeric(as.factor(std$group)),main = substitute(paste(italic('in vivo'), " SPME samples(selected peaks)")))
plotpca(std$data[corcluster$stdmassindex,],lv = as.numeric(as.factor(std$group)),main = substitute(paste(italic('in vivo'), " SPME samples(selected peaks by correlationship)")))

## ----globalcor-----------------------------------------------------------
std2 <- getstd(pmd,corcutoff = 0.9)
par(mfrow = c(1,3),mar = c(4,4,2,1)+0.1)
plotpca(std2$data,lv = as.numeric(as.factor(std2$group)),main = substitute(paste(italic('in vivo'), " SPME samples(all peaks)")))
plotpca(std$data[std$stdmassindex,],lv = as.numeric(as.factor(std$group)),main = substitute(paste(italic('in vivo'), " SPME samples(selected peaks)")))
plotpca(std2$data[std2$stdmassindex,],lv = as.numeric(as.factor(std2$group)),main = substitute(paste(italic('in vivo'), " SPME samples(selected peaks)")))

## ----sda-----------------------------------------------------------------
sda <- getsda(std, freqcutoff = 10)

## ----stdsda--------------------------------------------------------------
plotstdsda(sda)

## ----stdsdaidx-----------------------------------------------------------
par(mfrow = c(2,3),mar = c(4,4,2,1)+0.1)
plotstdsda(sda,sda$sda$diff2 == 0)
plotstdsda(sda,sda$sda$diff2 == 13.98)
plotstdsda(sda,sda$sda$diff2 == 15.99)
plotstdsda(sda,sda$sda$diff2 == 14.02)
plotstdsda(sda,sda$sda$diff2 == 28.03)
plotstdsda(sda,sda$sda$diff2 == 58.04)

## ----all,eval=F----------------------------------------------------------
#  sdaall <- getsda(spmeinvivo)
#  par(mfrow = c(2,3),mar = c(4,4,2,1)+0.1)
#  plotstdsda(sdaall,sdaall$sda$diff2 == 0)
#  plotstdsda(sdaall,sdaall$sda$diff2 == 13.98)
#  plotstdsda(sdaall,sdaall$sda$diff2 == 15.99)
#  plotstdsda(sdaall,sdaall$sda$diff2 == 14.02)
#  plotstdsda(sdaall,sdaall$sda$diff2 == 28.03)
#  plotstdsda(sdaall,sdaall$sda$diff2 == 58.04)

## ----rda-----------------------------------------------------------------
sda <- getrda(spmeinvivo$mz[std$stdmassindex])

## ----tarnet--------------------------------------------------------------
library(igraph)
# check metabolites of C18H39NO
chain <- getchain(spmeinvivo,diff = c(2.02,14.02,15.99,58.04,13.98),mass = 286.3101)
# show as network
net <- graph_from_data_frame(chain$sdac,directed = F)
plot(net,vertex.label=NA,vertex.size = as.numeric((names(V(net))))/25,edge.width = E(net)$diff2/20+0.00001)

## ----net-----------------------------------------------------------------
sda <- getsda(std, freqcutoff = 10)
df <- sda$sda
net <- graph_from_data_frame(df,directed = F)
plot(net,vertex.label=NA,vertex.size = as.numeric((names(V(net))))/25,edge.width = E(net)$diff2/20+0.00001)
# Check the degree of the nodes
# Show the degree distribution of the vertices
deg <- degree(net, mode="all")
degree_distribution(net)
plot(net, vertex.size=deg/2,vertex.label=NA,vertex.size = as.numeric((names(V(net))))/25, edge.width = E(net)$diff2/20+0.00001)
# network community structure detection
ceb <- cluster_edge_betweenness(net,weights = abs(E(net)$cor), directed = F) 
plot(ceb, net,vertex.label=NA,) 

## ----wrap----------------------------------------------------------------
result <- globalstd(spmeinvivo)


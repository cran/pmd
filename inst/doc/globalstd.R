## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----demodata-----------------------------------------------------------------
library(pmd)
data("spmeinvivo")
str(spmeinvivo)

## ----rtg----------------------------------------------------------------------
pmd <- getpaired(spmeinvivo, rtcutoff = 10, ng = 10)
plotrtg(pmd)

## ----pmd----------------------------------------------------------------------
plotpaired(pmd)

## ----pmdindex-----------------------------------------------------------------
# show the unique PMD found by getpaired function
for(i in 1:length(unique(pmd$paired$diff2))){
        diff <- unique(pmd$paired$diff2)[i]
        index <- pmd$paired$diff2 == diff
        plotpaired(pmd,index)
}

## ----std----------------------------------------------------------------------
std <- getstd(pmd)

## ----stdplot------------------------------------------------------------------
plotstd(std)

## ----stdrtplot----------------------------------------------------------------
par(mfrow = c(2,3))
plotstdrt(std,rtcluster = 23,main = 'Retention time group 23')
plotstdrt(std,rtcluster = 9,main = 'Retention time group 9')
plotstdrt(std,rtcluster = 18,main = 'Retention time group 18')
plotstdrt(std,rtcluster = 67,main = 'Retention time group 67')
plotstdrt(std,rtcluster = 49,main = 'Retention time group 49')
plotstdrt(std,rtcluster = 6,main = 'Retention time group 6')

## ----gettarget----------------------------------------------------------------
# you need retention time for independent peaks
index <- gettarget(std$rt[std$stdmassindex])
# output the ions for each injection
table(index)
# show the ions for the first injection
std$mz[index==1]
std$rt[index==1]

## ----pca----------------------------------------------------------------------
library(enviGCMS)
par(mfrow = c(1,2),mar = c(4,4,2,1)+0.1)
plotpca(std$data,lv = as.numeric(as.factor(std$group)),main = substitute(paste(italic('in vivo'), " SPME samples(all peaks)")))
plotpca(std$data[std$stdmassindex,],lv = as.numeric(as.factor(std$group)),main = substitute(paste(italic('in vivo'), " SPME samples(selected peaks)")))

## ----comp---------------------------------------------------------------------
stdcluster <- getcluster(std)
# extract pseudospectra for std peak 71
idx <- unique(stdcluster$cluster$largei[stdcluster$cluster$i==71])
plot(stdcluster$cluster$mz[stdcluster$cluster$largei==idx],stdcluster$cluster$ins[stdcluster$cluster$largei==idx],type = 'h',xlab = 'm/z',ylab = 'intensity',main = 'pseudospectra for GlobalStd peak 71')
# export peaks with the highest intensities in each GlobalStd peaks groups.
data <- stdcluster$data[stdcluster$stdmassindex2,]

## ----corcomp------------------------------------------------------------------
corcluster <- getcorcluster(spmeinvivo)
par(mfrow = c(1,3),mar = c(4,4,2,1)+0.1)
plotpca(std$data,lv = as.numeric(as.factor(std$group)),main = substitute(paste(italic('in vivo'), " SPME samples(all peaks)")))
plotpca(std$data[std$stdmassindex,],lv = as.numeric(as.factor(std$group)),main = substitute(paste(italic('in vivo'), " SPME samples(selected peaks)")))
plotpca(std$data[corcluster$stdmassindex,],lv = as.numeric(as.factor(std$group)),main = substitute(paste(italic('in vivo'), " SPME samples(selected peaks by correlationship)")))

## ----globalcor----------------------------------------------------------------
std2 <- getstd(pmd,corcutoff = 0.9)
par(mfrow = c(1,3),mar = c(4,4,2,1)+0.1)
plotpca(std2$data,lv = as.numeric(as.factor(std2$group)),main = substitute(paste(italic('in vivo'), " SPME samples(all peaks)")))
plotpca(std$data[std$stdmassindex,],lv = as.numeric(as.factor(std$group)),main = substitute(paste(italic('in vivo'), " SPME samples(selected peaks)")))
plotpca(std2$data[std2$stdmassindex,],lv = as.numeric(as.factor(std2$group)),main = substitute(paste(italic('in vivo'), " SPME samples(selected peaks)")))

## ----sda----------------------------------------------------------------------
sda <- getsda(std)

## ----stdsda-------------------------------------------------------------------
plotstdsda(sda)

## ----stdsdaidx----------------------------------------------------------------
par(mfrow = c(1,3),mar = c(4,4,2,1)+0.1)
plotstdsda(sda,sda$sda$diff2 == 2.02)
plotstdsda(sda,sda$sda$diff2 == 28.03)
plotstdsda(sda,sda$sda$diff2 == 58.04)

## ----all----------------------------------------------------------------------
sdaall <- getsda(spmeinvivo)
par(mfrow = c(1,3),mar = c(4,4,2,1)+0.1)
plotstdsda(sdaall,sdaall$sda$diff2 == 2.02)
plotstdsda(sdaall,sdaall$sda$diff2 == 28.03)
plotstdsda(sdaall,sdaall$sda$diff2 == 58.04)

## ----rda----------------------------------------------------------------------
sda <- getrda(spmeinvivo$mz[std$stdmassindex])

## ----tarnet-------------------------------------------------------------------
library(igraph)
# check metabolites of C18H39NO
chain <- getchain(spmeinvivo,diff = c(2.02,14.02,15.99,58.04,13.98),mass = 286.3101,digits = 2,corcutoff = 0)
# show as network
net <- graph_from_data_frame(chain$sdac,directed = F)
pal <- (grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11,"RdYlBu")
                )))(5)
plot(net,vertex.label=round(as.numeric(V(net)$name),2),vertex.size =5,edge.width = 5,edge.color = pal[as.numeric(as.factor(E(net)$diff2))],vertex.label.dist=1,vertex.color=ifelse(round(as.numeric(V(net)$name),4) %in% 286.3101,'red','black'),main = 'PMD network')
legend("topright",bty = "n",
       legend=unique(E(net)$diff2),
       fill=unique(pal[as.numeric(as.factor(E(net)$diff2))]), border=NA,horiz = F)

## ----net----------------------------------------------------------------------
sda <- getsda(std)
df <- sda$sda
net <- graph_from_data_frame(df,directed = F)
pal <- (grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11,"RdYlBu")
                )))(length(unique(E(net)$diff2)))
plot(net,vertex.label=round(as.numeric(V(net)$name)),vertex.size = 7,edge.width = 5,edge.color = pal[as.numeric(as.factor(E(net)$diff2))],main = 'PMD network')
legend("topright",bty = "n",
       legend=unique(E(net)$diff2),
       fill=unique(pal[as.numeric(as.factor(E(net)$diff2))]), border=NA,horiz = F)
# Check the degree of the nodes
# Show the degree distribution of the vertices
deg <- degree(net, mode="all")
degree_distribution(net)
plot(net, vertex.size=deg/2,vertex.label=NA,vertex.size = 7, edge.width = 5)
# network community structure detection
ceb <- cluster_edge_betweenness(net,weights = abs(E(net)$cor), directed = F) 
plot(ceb, net,vertex.label=NA,) 

## ----source-------------------------------------------------------------------
median(deg)
endogenous <- names(deg)[deg>median(deg)]
exogenous <- names(deg)[deg<=median(deg)]

## ----wrap---------------------------------------------------------------------
result <- globalstd(spmeinvivo,ng=10)

## -----------------------------------------------------------------------------
# all reaction
data("omics")
head(omics)
# kegg reaction
data("keggrall")
head(keggrall)
# literature reaction for mass spectrometry
data("sda")
head(sda)

## -----------------------------------------------------------------------------
data("hmdb")
head(hmdb)


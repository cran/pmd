## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(pmd)
data("spmeinvivo")

## ----tarnet-------------------------------------------------------------------
library(igraph)
# check metabolites of C18H39NO
# Use common PMDs for biological reactions
chain <- getchain(spmeinvivo,diff = c(2.02,14.02,15.99,58.04,13.98),mass = 286.3101,digits = 2,corcutoff = 0)
# show as network
net <- graph_from_data_frame(chain$sdac,directed = F)
pal <- grDevices::rainbow(5)
plot(net,vertex.label=round(as.numeric(V(net)$name),2),vertex.size =5,edge.width = 3,edge.color = pal[as.numeric(as.factor(E(net)$diff2))],vertex.label.dist=1,vertex.color=ifelse(round(as.numeric(V(net)$name),4) %in% 286.3101,'red','black'), main = 'PMD network')
legend("topright",bty = "n",
       legend=unique(E(net)$diff2),
       fill=unique(pal[as.numeric(as.factor(E(net)$diff2))]), border=NA,horiz = F)
# Consider the correlation coefficient cutoff
chain <- getchain(spmeinvivo,diff = c(2.02,14.02,15.99,58.04,13.98),mass = 286.3101,digits = 2,corcutoff = 0.6)
# show as network
net <- graph_from_data_frame(chain$sdac,directed = F)
pal <- grDevices::rainbow(5)
plot(net,vertex.label=round(as.numeric(V(net)$name),2),vertex.size =5,edge.width = 3,edge.color = pal[as.numeric(as.factor(E(net)$diff2))],vertex.label.dist=1,vertex.color=ifelse(round(as.numeric(V(net)$name),4) %in% 286.3101,'red','black'), main = 'PMD network')
legend("topright",bty = "n",
       legend=unique(E(net)$diff2),
       fill=unique(pal[as.numeric(as.factor(E(net)$diff2))]), border=NA,horiz = F)

## ----net----------------------------------------------------------------------
std <- globalstd(spmeinvivo,sda = F)
sda <- getsda(std,freqcutoff = 12)
df <- sda$sda
net <- graph_from_data_frame(df,directed = F)
pal <- grDevices::rainbow(length(unique(E(net)$diff2)))
plot(net,vertex.label=NA,vertex.size = 5,edge.width = 3,edge.color = pal[as.numeric(as.factor(E(net)$diff2))],main = 'PMD network')
legend("topright",bty = "n",
       legend=unique(E(net)$diff2),
       fill=unique(pal[as.numeric(as.factor(E(net)$diff2))]), border=NA,horiz = F)

## ----nwa----------------------------------------------------------------------
# network community structure detection
ceb <- cluster_edge_betweenness(net,weights = abs(E(net)$cor), directed = F) 
plot(ceb, net,vertex.label=NA,vertex.size = 5,edge.width = 3,) 
# output membership
head(cbind(ceb$membership,ceb$names))

## ----cda----------------------------------------------------------------------
cbp <- enviGCMS::getfilter(std,rowindex = std$stdmassindex)
cda <- getcda(cbp)
df <- cda$cda
# filter based on retention time differences larger than 2 mins
df <- df[df$diffrt>120,]
netc <- graph_from_data_frame(df,directed = F)
plot(netc,vertex.label=NA,vertex.size = 5,edge.width = 3,main = 'Correlation network')

## ----source-------------------------------------------------------------------
deg <- degree(net, mode = 'all')
median(deg)
endogenous <- names(deg)[deg>median(deg)]
exogenous <- names(deg)[deg<=median(deg)]

## -----------------------------------------------------------------------------
pmd <- getreact(spmeinvivo,pmd=15.99,ratiocv = 30)
data15.99Da <- apply(pmd$data,2,sum)
# show the quantative reaction level across samples
data15.99Da

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

## -----------------------------------------------------------------------------
plotcn('C6H12O6','Glucose',c(2.016,14.016,15.995))


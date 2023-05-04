#!/usr/bin/env Rscript

library(igraph)
library(stringr)
library(dplyr)
#library(readr)
library(svglite)

### great resource: https://users.dimi.uniud.it/~massimo.franceschet/R/communities.html
### the original documentation for the cluster_walktrap algorithm https://igraph.org/r/doc/cluster_walktrap.html
### the idea is that in a graph (here, nodes are connected by their angstrom distance edges) short walks along those 
### edges will be able to discriminate community memberships, since short walks from starting point X should be
### staying in the same community (as point X is in)

args = commandArgs(trailingOnly=TRUE)
################################# C-beta Version ################################# ################################# 
filus <- args[1]
namus <- args[2]
outfilus <- args[3]
svgfilus <- args[4]
################ 
newdat <- read.csv(filus,row.names=1, quote="", stringsAsFactors=FALSE, check.names = FALSE)
maintilus <- c(namus, "Functional residues\nwith Cb distances below 8A")
#View(newdat)
rownames(newdat) <- colnames(newdat)
data <- as.matrix(as.data.frame(lapply(newdat, as.numeric),check.names = FALSE))
netw <- graph_from_adjacency_matrix(data)
#### removing unconnected edges
Isolated = which(degree(netw)==0)     
G2 = delete.vertices(netw, Isolated)
coords = layout_with_fr(G2)
#### plotting
c1 =cluster_walktrap(G2)
svglite(svgfilus, width = 10, height = 10)
plot(c1,
     G2,
     layout=coords,
     main=paste(maintilus, collapse = " "),
     vertex.size = 1,
     vertex.label.cex=0.7,
     vertex.label.color="black",
     vertex.frame.color="transparent",
     edge.arrow.size=0,
     edge.width=0.05
)

dev.off()

max(membership(c1))
for (si in sizes(c1)){
  #cat(si)
}


fileConn <- file(outfilus,"a")
for (i in 1:max(membership(c1))){
  clustersize <- nth(sizes(c1),i)
  if (clustersize == 1){
    residue <- c1[[c(i,1)]]
    #cat(namus,"\t",residue,"\n")
  }
  else {
    for (a in 1:clustersize){
      residue <- c1[[c(i,a)]]      
      txt <- paste(namus,residue,i,sep="\t")
      #cat(namus,"\t",residue,"\t",i,"\n")
      #cat(txt)
      cat(txt, file = fileConn,sep="\n",append = TRUE)
    }
  }
}

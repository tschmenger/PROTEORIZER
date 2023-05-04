#!/usr/bin/env Rscript

library(stringr)
library(dplyr)
#library(readr)
library(svglite)
library(tidyr)

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

###############################
### hierarchical clustering
### preparing the data by introducing NA for every empty cell or cells with []

# hieradata <- lapply(newdat, function(x) as.numeric(as.character(x)))
# hieradata <- as.data.frame(hieradata)
# ### converting NAs to a big number, 100, to make them irrelevant
# hieradata[is.na(hieradata)] =100
# hieradata

newdat[newdat == "[]"] <- 99999
newdat[newdat == ""] <- 99999
dataset <- as.data.frame(Filter(function(x)(length(unique(x))>1), newdat))
dataset <- dataset[apply(dataset,
                         MARGIN = 1, # rowwise
                         FUN = function(x) length(unique(x)) > 1), ]

### calculating a distance matrix
# dist_hieradata <- dist(as.data.frame(hieradata))
distances <- dist(dataset,method = "euclidean")

### performing the hierarchical clustering
# hieraclusters <- hclust(dist_hieradata)
# plot(hieraclusters)
hieraclusters <- hclust(distances)

### automatically determine the number of clusters we got
# library(dynamicTreeCut)
# betterclusternumber <- cutreeDynamicTree(hieraclusters, minModuleSize = 3, deepSplit = FALSE)
# number_to_cut <- length(which(betterclusternumber != 0))

### cutting the clustered tree above to give us "betterclusternumber" clusters
# clusterCut <- cutree(hieraclusters, k = number_to_cut)
height_to_cut <- max(hieraclusters$height) - (0.25 * (max(hieraclusters$height) - min(hieraclusters$height)))
clusterCut <- cutree(hieraclusters, h = height_to_cut)



fileConn <- file(outfilus,"a")

for ( i in 1:length(colnames(dataset))){### getting rid of clusters of size < 3
  if (length(which(clusterCut == clusterCut[i]))>2){
  txt <- paste(namus, colnames(dataset[i]),clusterCut[i], sep="\t")
  cat(txt, file = fileConn,sep="\n",append = TRUE)
}}


library(ggplot2)
library(cowplot)

positions <- colnames(dataset)
hClust <- clusterCut

results <- data.frame(as.numeric(positions),
                      as.numeric(hClust))

### Hierarchical Clustering, Plotting
svglite(svgfilus, width = 10, height = 10)

hier <- ggplot(results,aes(x = as.integer(as.character(positions)), y = hClust))+
  geom_point(aes(color=as.factor(hClust)))+
  theme_classic()+
  ggtitle(paste(maintilus, collapse = " "))+ # for the main title 
  xlab(paste(namus,"Residue",sep=" "))+ # for the x axis label
  ylab("Cluster")+ # for the y axis label
  theme(axis.text.x=element_text(angle = 90, size = 8))+
  labs(color = "Clusters")
hier

dev.off()



library(ape)
library(adegenet)
library(knitr)
library(igraph)
library("RColorBrewer")


sym <- function(M) {
  M[lower.tri(M)] = t(M)[lower.tri(M)]
  M
}

#args = commandArgs(trailingOnly=TRUE)
#distmatF <- args[1]


distmatF <- "Thies_all_manual.PASS.Cls.miss0.5.LMRG.HAP.INDEL.recode.vcf.dist.tab.txt" 
mat <- read.table(distmatF,sep="\t")
D <- as.dist(sym(mat))

#make into clusters
clust <- gengraph(D,ngrp=3)
#plot(clust$g, main="gengraph clusters")
distmat <- as.matrix(D)
names <- colnames(distmat)
years <- as.numeric(gsub("^.*\\.","",names))
years <- years+2000
colls <- as.Date(paste("1","jan",years,sep=""),"%d%b%Y")
names(colls)<-names
  
weeks <- as.integer(difftime(coll, min(coll), unit="weeks"))

for (cID in unique(clust$clust$membership)) {
  name <- names[clust$clust$membership==cID]
  year <- years[clust$clust$membership==cID]
  coll <- colls[name]
  dist <- distmat[name,name] 
  res <- seqTrack(dist, x.names=name, x.dates=coll)
  
  span <- max(years)-min(years)+1
  cols <- rev(brewer.pal(span, "RdBu") )
  
  ts=1 #textsize
  ig <- as.igraph(res)
  V(ig)$name <- name
  V(ig)$year <- year
  V(ig)$color <- cols[year-min(year)+1]
#  V(ig)$label.cex <- ts
  write.graph(ig,paste(distmatF,cID,"dot",sep="."),format="dot")
  write.graph(ig,paste(distmatF,cID,"gml",sep="."),format="graphml")

  plot(ig,main="all vars",vertex.size=25)
}
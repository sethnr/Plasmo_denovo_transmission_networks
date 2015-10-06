library(ggplot2)
library(reshape2)


fileroot <- commandArgs(TRUE)[1]

pos <- read.table(paste(fileroot,"012","pos",sep="."))
indv <- read.table(paste(fileroot,"012","indv",sep="."))[,1]
genos <- t(read.table(paste(fileroot,"012",sep=".")))
genos[genos==-1] <- NA
colnames(genos)<-indv
genotab <- cbind(pos,genos[2:dim(genos)[[1]],])
colnames(genotab)[1:2]<-c("chr","pos")
genotab.m <- melt(genotab,id.vars =c("chr","pos"),variable.name ="indv")
genotab.m$value <- as.factor(genotab.m$value)
#sort levels
genotab.m$indv <- factor(genotab.m$indv,levels=sort(levels(genotab.m$indv)))

ggplot(genotab.m,aes(x=pos,y=indv,colour=value,group=indv)) + geom_point() + facet_grid(chr ~ .)

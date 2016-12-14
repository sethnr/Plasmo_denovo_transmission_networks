library(ggplot2)
library(ape)
library(stringdist)
library(reshape2)

fileroot<-"Thies_all_manual.PASS.Cls.miss0.5.LMRG.HAP.alleles"
fileroot <- commandArgs(TRUE)[1]

tab=paste(fileroot,'tab',sep='.')
outpng=paste(fileroot,'png',sep='.')
outnex=paste(fileroot,'nexus',sep='.')
outtab=paste(fileroot,'dist.tab.txt',sep='.')

#read in ped file and get hamming distances:
alleleTab <- read.table(tab,colClasses="character",header=T,na.strings = c("."))

genos <- t(data.matrix(alleleTab[6:dim(alleleTab)[2]]))
inds <- row.names(genos)

genosDat <- as.phyDat(genos, type="USER", levels = c(0:max(genos,na.rm=T)))

distmat <- dist.hamming(genosDat)
treeNJ <- NJ(distmat)


genos[1:2,genos[1,]!=genos[2,]]

#
distmat2 <- matrix(nrow=length(inds),ncol=length(inds))
colnames(distmat2) = inds
rownames(distmat2) = inds

write("calculating distance matrix",stderr())
 for (i in rownames(genos)){
   for (j in rownames(genos)){
      notnull <- !(is.na(genos[i,]) | is.na(genos[j,]))
      
      distmat2[i,j] <- sum(genos[i,notnull] != genos[j,notnull])
  }
}

#get parsinomy score
parsimony(treeNJ,genosDat)

#get parsinomy tree
#treePars <- optim.parsimony(treeNJ,genosDat)
treeRatchet <- pratchet(genosDat)
#re-add distances
treeRatchet <- nnls.phylo(treeRatchet,distmat)
treeRatchet <- midpoint(treeRatchet)

#parsimony(c(treePars,treeRatchet),genosDat)
#plot(treePars,"unrooted",main="parsimony: NNI swaps")
png(filename=outpng,type="cairo-png")
plot(treeRatchet,show.node.label = T,show.tip.label = T,
     main=paste("parsimony ratchet, ",fileroot),type="phylogram",cex = 1.1)
dev.off()


write.nexus(treeRatchet,file=outnex)

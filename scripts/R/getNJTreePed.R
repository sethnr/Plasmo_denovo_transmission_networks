library(ggplot2)
library(ape)
library(stringdist)
library(reshape2)


fileroot <- commandArgs(TRUE)[1]

ped=paste(fileroot,'ped',sep='.')
outpng=paste(fileroot,'png',sep='.')
outnex=paste(fileroot,'nexus',sep='.')
outtab=paste(fileroot,'dist.tab.txt',sep='.')

#read in ped file and get hamming distances:
genos <- read.table(ped,colClasses="character")
inds <- genos[,1]
genos <- genos[,seq(7,dim(genos)[[2]],2)]
rownames(genos)=inds

distmat = matrix(nrow=length(inds),ncol=length(inds))
colnames(distmat) = inds
rownames(distmat) = inds

write("calculating distance matrix",stderr())
for (i in rownames(genos)){
    for (j in rownames(genos)){
    filled = intersect(which(genos[i,] != 0), which(genos[j,] !=0))
    write(paste("calculating",i,"v",j),stderr())
    write(length(filled),stderr())
    
#    distmat[i,j] = sum(genos[i,]!=genos[j,])
    distmat[i,j] = sum(genos[i,filled]!=genos[j,filled])
    }
}

write.table(distmat,stderr())
            
#make NJ trees and save as Nexus
njtree <- nj(as.dist(distmat))

Sys.setenv("DISPLAY"=":0.0")

png(filename=outpng,type="cairo-png")
plot(njtree,show.node.label = T,show.tip.label = T,
     main=paste("NJ tree, ",fileroot),type="fan",cex = 1.1) 
dev.off()

write.nexus(njtree,file=outnex)

distmat[lower.tri(distmat)] <- ''
write.table(distmat,sep="\t",file=outtab,quote=F)

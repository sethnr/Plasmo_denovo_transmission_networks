library(ggplot2)
library(ape)
library(stringdist)
library(reshape2)


fileroot <- commandArgs(TRUE)[1]
subSample=1  #proportion to use for subsample

boot <- Sys.getenv('SGE_TASK_ID')

ped=paste(fileroot,'ped',sep='.')
outpng=paste(fileroot,boot,'png',sep='.')
outnex=paste(fileroot,boot,'nexus',sep='.')
outtab=paste(fileroot,boot,'dist.tab.txt',sep='.')

#read in ped file and get hamming distances:
genos <- read.table(ped,colClasses="character")


inds <- genos[,1]
genos <- genos[,seq(7,dim(genos)[[2]],2)]
noGenos = dim(genos)[[2]]


#sample <- sample(1:noGenos,round(noGenos*subSample))
sample <- sample(1:noGenos,replace=T)

rownames(genos)=inds
genos <- genos[,sample]



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
    if(length(filled) > 0) {
        distmat[i,j] = sum(genos[i,filled]!=genos[j,filled])
      }else {
    	  distmat[i,j]=-1
      }
    }
}

write.table(distmat,stderr())
            
#make NJ trees and save as Nexus
njtree <- nj(as.dist(distmat))

#Sys.setenv("DISPLAY"=":0.0")

#png(filename=outpng,type="cairo-png")
#plot(njtree,show.node.label = T,show.tip.label = T,
#     main=paste("NJ tree, ",fileroot),type="fan",cex = 1.1) 
#dev.off()

write.nexus(njtree,file=outnex)

#distmat[lower.tri(distmat)] <- ''
#write.table(distmat,sep="\t",file=outtab,quote=F)

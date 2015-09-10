#library(PopGenome)
library(ggplot2)
library(ape)
library(stringdist)
library(reshape2)

#read in ped file and get hamming distances:
genos <- read.table("GATK_calls_daniels_Thies.ped",colClasses="character")
inds <- genos[,1]
genos <- genos[,seq(7,dim(genos)[[2]],2)]
rownames(genos)=inds

distDisco = matrix(nrow=length(inds),ncol=length(inds))
colnames(distDisco) = inds
rownames(distDisco) = inds

for (i in rownames(genos)){
    for (j in rownames(genos)){
#    write(i,stderr())
    distDisco[i,j] = sum(genos[i,]!=genos[j,])
    }
}
            

#read in ped file and get hamming distances:
genos <- read.table("daniels_2015.tab.txt",colClasses="character",sep="\t",skip=1)
inds <- genos[,1]
genos <- genos[,4:dim(genos)[[2]]]
rownames(genos)=inds

distDan = matrix(nrow=length(inds),ncol=length(inds))
colnames(distDan) = inds
rownames(distDan) = inds

#distdan<-stringdistmatrix(genos,method="hamming")

n=0
for (i in rownames(genos)){
    n = n+1
    write(paste(n,i),stderr())
    for (j in rownames(genos)){
    distDan[i,j] = sum(genos[i,]!=genos[j,])
    }
}
            


#make NJ trees and save as Nexus
danielsNJ <- nj(as.dist(distDan))
plot(danielsNJ,show.node.label = T,show.tip.label = T,
     main="NJ tree, barcode SNPs, Daniels '15",type="radial",cex = 1.1) 

write.nexus(danielsNJ,file="daniels.thies.nexus")

discoNJ <- nj(as.dist(distDisco))
plot(discoNJ,show.node.label = T,show.tip.label = T,
     main="NJ tree, barcode SNPs, Discovar samples",type="cladogram",cex = 1.1) 

write.nexus(discoNJ,file="disco.thies.nexus")

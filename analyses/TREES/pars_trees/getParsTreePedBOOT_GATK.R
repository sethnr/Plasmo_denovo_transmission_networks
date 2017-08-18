#library(ggplot2)
library(phangorn)
library(ape)

sym <- function(M) {
  M[lower.tri(M)] = t(M)[lower.tri(M)]
  M
}



#fileroot<-"Thies_all_manual.PASS.Cls.miss0.5.LMRG.HAP.DISCORDS.miss-1.alleles"
#outfolder <- "boot_discords"

fileroot<-"thies_300100_haplo.CALLBOTH.RENAME.vcf.gz"
outfolder <- "boot_gatk_both"

#fileroot <- commandArgs(TRUE)[1]


tab=paste(fileroot,'tab',sep='.')
#boot <- Sys.getenv('SGE_TASK_ID')


for (boot in c(1:100)) {
    
  
  outnex=paste(".",outfolder, paste(fileroot,boot,'pars.nexus',sep='.'),sep="/")
  outtab=paste(".",outfolder, paste(fileroot,boot,'dist.tab.txt',sep='.')         ,sep="/")
  
  
  #read in ped file and get hamming distances:
  alleleTab <- read.table(tab,colClasses="character",header=T,na.strings = c("."))
  
  genos <- t(data.matrix(alleleTab[6:dim(alleleTab)[2]]))
  inds <- row.names(genos)
  
  
  noGenos = dim(genos)[[2]]
  sample <- sample(1:noGenos,replace=T)
  
  #bootstrap to full size of genotype matrix:
  genos <- genos[,sample]
  
  genosDat <- as.phyDat(genos, type="USER", levels = c(0:max(genos,na.rm=T)))
  
  distmat <- dist.hamming(genosDat)
  #treeNJ <- NJ(distmat)
  
  #
  #distmat2 <- matrix(nrow=length(inds),ncol=length(inds))
  #colnames(distmat2) = inds
  #rownames(distmat2) = inds
  #write("calculating distance matrix",stderr())
  #for (i in rownames(genos)){
  #   for (j in rownames(genos)){
  #      notnull <- !(is.na(genos[i,]) | is.na(genos[j,]))
  #      
  #      distmat2[i,j] <- sum(genos[i,notnull] != genos[j,notnull])
  #  }
  #}
  
  
  #get parsinomy score
  #parsimony(treeNJ,genosDat)
  
  #get parsinomy tree
  #treePars <- optim.parsimony(treeNJ,genosDat)
  write("getting parsimony tree",stderr())
  treeRatchet <- pratchet(genosDat)
  write(str(treeRatchet))
  #re-add distances
#  write(dim(as.matrix(distmat)),stderr())
#  write.table(round(as.matrix(distmat),2),stderr())
  #write(treeRatchet$tip.label,stderr())
  
  
  write("re-adding distances",stderr())
  treeRatchet <- nnls.phylo(treeRatchet,as.matrix(distmat),rooted=T)
  #treeRatchet <- nnls.phylo(treeRatchet,distmat)
  
  #write(treeRatchet,stderr())
  
  #DO NOT MIDPOINT ROOT FOR BOOTSTRAPPING
  #write("midpoint root",stderr())
#  treeRatchet <- midpoint(treeRatchet)
  
  
  #write("printing tree",stderr())
  #png(filename=outpng,type="cairo-png")
#  plot(treeRatchet,show.node.label = T,show.tip.label = T,
#       main=paste("parsimony ratchet, ",fileroot),type="phylogram",cex = 1.1)
  #dev.off()
  write.nexus(treeRatchet,file=outnex)
}



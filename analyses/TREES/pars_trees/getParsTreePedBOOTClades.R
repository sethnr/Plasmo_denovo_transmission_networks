#library(ggplot2)

sym <- function(M) {
  M[lower.tri(M)] = t(M)[lower.tri(M)]
  M
}



fileroot<-"Thies_all_manual.PASS.Cls.miss0.5.LMRG.HAP.DISCORDS.miss-1.alleles"
outfolder<-"boot_discords_clades"
#fileroot <- commandArgs(TRUE)[1]

tab=paste(fileroot,'tab',sep='.')
boot <- Sys.getenv('SGE_TASK_ID')

c1 <- c("Th086.07","Th106.09","Th106.11","Th117.11","Th132.11","Th134.11","Th162.12","Th196.12","Th230.12","Th074.13")
c2 <- c("Th166.12","Th092.13","Th211.13","Th245.13","Th246.13" )
c3 <- c("Th068.12","Th061.13","Th095.13")
c1o <- "Th245.13"
c2o <- "Th086.07"
c3o <- "Th086.07"
c1dis <- apply(alleleTab[,c1],1,function(x){length(unique(na.omit(x)))})>1
c2dis <- apply(alleleTab[,c2],1,function(x){length(unique(na.omit(x)))})>1
c3dis <- apply(alleleTab[,c3],1,function(x){length(unique(na.omit(x)))})>1

clades <- list(list("name"="c29","clade"=c1,"vars"=c1dis,"out"=c1o),
            list("name"="c26","clade"=c2,"vars"=c2dis,"out"=c2o),
            list("name"="c23","clade"=c3,"vars"=c3dis,"out"=c3o))

for (boot in c(1:100)) {
    
  
  
  for (clade in clades) {
#    ped=paste(fileroot,'ped',sep='.')
#    outpng=paste(fileroot,boot,'png',sep='.')
    
#    outnex=paste("./boot_discords_clades", paste(fileroot,clade$name,boot,'pars.nexus',sep='.'),sep="/")
#    outtab=paste("./boot_discords_clades", paste(fileroot,boot,'dist.tab.txt',sep='.')         ,sep="/")
    outnex=paste(".",outfolder, paste(fileroot,clade$name,boot,'pars.nexus',sep='.'),sep="/")
    outtab=paste(".",outfolder, paste(fileroot,boot,'dist.tab.txt',sep='.')         ,sep="/")
    
  #read in ped file and get hamming distances:
    alleleTab <- read.table(tab,colClasses="character",header=T,na.strings = c("."))
    genos <- t(alleleTab[clade$vars,c(clade$clade,clade$out)])
    
    
    #genos <- t(data.matrix(alleleTab[6:dim(alleleTab)[2]]))
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
    #write(dim(as.matrix(distmat)),stderr())
    #write.table(round(as.matrix(distmat),2),stderr())
    #write(treeRatchet$tip.label,stderr())
    
    
    write("re-adding distances",stderr())
    treeRatchet <- nnls.phylo(treeRatchet,as.matrix(distmat),rooted=T)
    #treeRatchet <- nnls.phylo(treeRatchet,distmat)
    
    
    write("outgroup root",stderr())
    treeRatchet <- root(treeRatchet,clade$out)
#     write("midpoint root",stderr())
#     treeRatchet <- midpoint(treeRatchet)
    
    #parsimony(c(treePars,treeRatchet),genosDat)
    #plot(treePars,"unrooted",main="parsimony: NNI swaps")
    
    
    #write("printing tree",stderr())
    #png(filename=outpng,type="cairo-png")
    #plot(treeRatchet,show.node.label = T,show.tip.label = T,
    #     main=paste("parsimony ratchet, ",fileroot),type="phylogram",cex = 1.1)
    #dev.off()
    write.nexus(treeRatchet,file=outnex)
  }
}



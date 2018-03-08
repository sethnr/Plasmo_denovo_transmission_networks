
```r
library(ape)
library(adegenet)
library(phangorn)
library(knitr)
#library(igraph)
library(RColorBrewer)
library(ggplot2)
library(reshape2)

require(gridExtra)


opts_chunk$set(fig.width=9, fig.height=9)
#opts_chunk$set(dev=c('png','postscript'))
opts_chunk$set(dev=c('png'))
```



```r
sym <- function(M) {
  M[lower.tri(M)] = t(M)[lower.tri(M)]
  M
}


makeDist <- function(distance_matrix_file, meta_file, ngroups=3) {
  mat <- read.table(distance_matrix_file,sep="\t")
  D <- as.dist(sym(mat))
  clust <- gengraph(D,ngrp=ngroups)
  names <- colnames(mat)
  mat <- as.matrix(mat)
  
  name1 <- names[clust$clust$membership==1]
  name2 <- names[clust$clust$membership==2]
  name3 <- names[clust$clust$membership==3]
  dist1 <- mat[name1,name1]
  dist2 <- mat[name2,name2]
  dist3 <- mat[name3,name3]
    
  list(dist1,dist2,dist3)
}
```



```r
getSplitSupports <- function(tree,genos) {
    
#  samps <- c(1:dim(genos)[[1]])
  samps <- tree$tip.label
  
  splitSupports <- data.frame(splits=character(),
                              null=numeric(),
                              irrelevant=numeric(),
                              pro=numeric(),
                              anti=numeric(),
                              stringsAsFactors = F)
  
  splits <- as.splits(tree)
  write(paste(c((length(samps)+1):length(splits)),sep=", "),stderr())
  for (si in c((length(samps)+1):length(splits))) {
    split <- samps[splits[[si]]]
    outs <- samps[!samps %in% split]
    #if (length(split)==1) {next}
    if (length(split)==length(samps)) {next}
    
  #  write(split,stderr())
    calcs <- apply(genos,2,function(x) {
      sAll <- na.omit(unique(x[split]))
      oAll <- na.omit(unique(x[outs]))
      comm <- intersect(sAll,oAll)
      sLen <- length(sAll)
      oLen <- length(oAll)
      cLen <- length(comm)
      if (sLen > 2) {sLen<-2}
      if (oLen > 2) {oLen<-2}
      if (cLen > 1) {cLen<-1}
      cat <- paste(c(sLen,cLen,oLen),sep="",collapse="")
      if (sLen==0 | oLen==0 | cat =="111") {
        return("null") #if null or monomorphic
      }
      if (cat %in% c("101","102","201","202")) {
        return("support")
      }
      if (cat %in% c("112","211")) {
        return("no support")
      }
      if (cat == "212") {
        return("against")
      }
      write(cat,stderr())
        return(cat)
      })
   
      #write(table(calcs),stderr()) 
      # write(paste(split,sep=",",collapse =","),stderr())
      # write(table(calcs)[c("null","no support","support","against")],stderr())
      nextI = dim(splitSupports)[[1]]+1
      splitSupports[nextI,"splits"] = paste(split,sep=",",collapse =",")
      counts <- table(calcs)[c("null","no support","support","against")]
      counts[is.na(counts)] <- 0
      splitSupports[nextI,c("null","irrelevant","pro","anti")] = counts
      
  }
  
  return(splitSupports)
}

getSplitPlot <- function(tree) {
    
#  samps <- c(1:dim(genos)[[1]])
  samps <- tree$tip.label
  
  splitPlot <- data.frame(splits=character(),
                          sample=character(),
                          insplit=numeric(),
                          stringsAsFactors = F)
  
  splits <- as.splits(tree)
  write(paste(c((length(samps)+1):length(splits)),sep=", "),stderr())
  for (si in c((length(samps)+1):length(splits))) {
    split <- samps[splits[[si]]]
    outs <- samps[!samps %in% split]
    #if (length(split)==1) {next}
    if (length(split)==length(samps)) {next}
    
    splitStr = paste(split,sep=",",collapse =",")
    
    write(paste(si,splitStr,length(splits)),stderr())
    for (s in split) {
      nextI = dim(splitPlot)[[1]]+1
      splitPlot[nextI,"splits"] = splitStr
      splitPlot[nextI,"sample"] = s
      splitPlot[nextI,"insplit"] = 1}
    for (s in outs) {
      nextI = dim(splitPlot)[[1]]+1
      splitPlot[nextI,"splits"] = splitStr
      splitPlot[nextI,"sample"] = s
      splitPlot[nextI,"insplit"] = 0}
    
      
  }
  
  treeOrder <- tree$tip.label[tree$edge[c(tree$edge[,2] <= length(tree$tip.label)),2]]
  splitPlot$sample <- factor(splitPlot$sample,levels=treeOrder,ordered=T)
  splitPlot$insplit <- as.logical(splitPlot$insplit)
  return(splitPlot)
}

blankTheme <- theme(axis.title=element_blank(),axis.text = element_blank())
blankThemeLab <- theme(axis.title=element_blank(),axis.text.y = element_blank(),
                    axis.text.x=element_text(angle=90,hjust=1))
```





```r
meta <- read.table("Thies_metadata_1701.txt",sep="\t",header=T)
colnames(meta)[1]<-"name"
meta <- meta[!is.na(meta$Age),]

discoDists <- read.table("thies_disco.FILT.m0.5.vcf.gz.dist.tab.txt",header=T,sep="\t")

meta <- subset(meta,name %in% names(discoDists))
rownames(meta) <- meta$name

hap100Dists <- read.table("thies_haplo100.FILT.m0.5.vcf.gz.dist.tab.txt",header=T,sep="\t")

hap250Dists <- read.table("thies_haplo250.FILT.m0.5.vcf.gz.dist.tab.txt",header=T,sep="\t")



samples <- read.table("../180224_core_callable_VQSLOD90/Thies.samples.txt",col.names=c("ID","lane","dataset","sample","bam"))
samples <- unique(samples[,c("ID","sample")])
snames <- samples[,c("ID")]; samples<- as.character(samples[,c("sample")]) ; names(samples) <- snames

#fix names for hap100 tree
names(hap100Dists) <- samples[rownames(hap100Dists)]
rownames(hap100Dists) <- samples[rownames(hap100Dists)]
```




```r
cl1 <- c("Th086.07", "Th106.09", "Th106.11", "Th117.11", "Th132.11", "Th134.11", "Th162.12", "Th196.12", "Th230.12", "Th074.13")
cl1yrs <- as.numeric(gsub(".*\\.","",cl1))
cl2 <- c("Th166.12", "Th092.13", "Th211.13", "Th245.13", "Th246.13")
cl3 <- c("Th068.12", "Th061.13", "Th095.13")
```


#MAKE DISCOVAR TREES

```r
discoAlleles <- read.table("thies_disco.FILT.m0.5.alleles.tab",colClasses="character",header=T,na.strings = c("."))
discoGenos <- t(data.matrix(discoAlleles[6:dim(discoAlleles)[2]]))
inds <- row.names(discoGenos)
genosDat <- as.phyDat(discoGenos, type="USER", levels = c(0:max(discoGenos,na.rm=T)))

# discoPTree <- pratchet(genosDat) # ,njtree1)
# discoPTree <- nnls.phylo(discoPTree,dist.hamming(genosDat))
# discoPTree <- midpoint(discoPTree)
discoPTree <- read.nexus("thies_disco.FILT.m0.5.vcf.gz.PARS.nexus")
discoPTree1 <- drop.tip(discoPTree,c(cl2,cl3))
discoPTree2 <- drop.tip(discoPTree,c(cl3,cl1))
discoPTree3 <- drop.tip(discoPTree,c(cl1,cl2))
plot(discoPTree)
```

![plot of chunk unnamed-chunk-6](figure/chrs-unnamed-chunk-6-1.png)

```r
plot(discoPTree1)
```

![plot of chunk unnamed-chunk-6](figure/chrs-unnamed-chunk-6-2.png)

```r
plot(discoPTree2)
```

![plot of chunk unnamed-chunk-6](figure/chrs-unnamed-chunk-6-3.png)

```r
plot(discoPTree3)
```

![plot of chunk unnamed-chunk-6](figure/chrs-unnamed-chunk-6-4.png)


#MAKE Hap250 TREES

```r
hap250Alleles <- read.table("thies_haplo250.FILT.m0.5.alleles.tab",colClasses="character",header=T,na.strings = c("."))
hap250Genos <- t(data.matrix(hap250Alleles[6:dim(hap250Alleles)[2]]))
inds <- row.names(hap250Genos)
genosDat <- as.phyDat(hap250Genos, type="USER", levels = c(0:max(hap250Genos,na.rm=T)))

# hap250PTree <- pratchet(genosDat) # ,njtree1)
# hap250PTree <- nnls.phylo(hap250PTree,dist.hamming(genosDat))
# hap250PTree <- midpoint(hap250PTree)
hap250PTree <- read.nexus("thies_haplo250.FILT.m0.5.vcf.gz.PARS.nexus")
hap250PTree1 <- drop.tip(hap250PTree,c(cl2,cl3))
hap250PTree2 <- drop.tip(hap250PTree,c(cl3,cl1))
hap250PTree3 <- drop.tip(hap250PTree,c(cl1,cl2))
plot(hap250PTree)
```

![plot of chunk unnamed-chunk-7](figure/chrs-unnamed-chunk-7-1.png)

```r
plot(hap250PTree1)
```

![plot of chunk unnamed-chunk-7](figure/chrs-unnamed-chunk-7-2.png)

```r
plot(hap250PTree2)
```

![plot of chunk unnamed-chunk-7](figure/chrs-unnamed-chunk-7-3.png)

```r
plot(hap250PTree3)
```

![plot of chunk unnamed-chunk-7](figure/chrs-unnamed-chunk-7-4.png)


#MAKE Hap100 TREES

```r
hap100Alleles <- read.table("thies_haplo100.FILT.m0.5.alleles.tab",colClasses="character",header=T,na.strings = c("."))
cnames <- gsub("SM.","SM-",colnames(hap100Alleles))
cnames[cnames %in% names(samples)] <- samples[cnames[cnames %in% names(samples)]]
colnames(hap100Alleles) <- cnames


hap100Genos <- t(data.matrix(hap100Alleles[6:dim(hap100Alleles)[2]]))
inds <- row.names(hap100Genos)
genosDat <- as.phyDat(hap100Genos, type="USER", levels = c(0:max(hap100Genos,na.rm=T)))

# hap100PTree <- pratchet(genosDat) # ,njtree1)
# hap100PTree <- nnls.phylo(hap100PTree,dist.hamming(genosDat))
# hap100PTree <- midpoint(hap100PTree)


hap100PTree <- read.nexus("thies_haplo100.FILT.m0.5.vcf.gz.PARS.nexus")
hap100PTree1 <- drop.tip(hap100PTree,c(cl2,cl3))
hap100PTree2 <- drop.tip(hap100PTree,c(cl3,cl1))
hap100PTree3 <- drop.tip(hap100PTree,c(cl1,cl2))
plot(hap100PTree)
```

![plot of chunk unnamed-chunk-8](figure/chrs-unnamed-chunk-8-1.png)

```r
plot(hap100PTree1)
```

![plot of chunk unnamed-chunk-8](figure/chrs-unnamed-chunk-8-2.png)

```r
plot(hap100PTree2)
```

![plot of chunk unnamed-chunk-8](figure/chrs-unnamed-chunk-8-3.png)

```r
plot(hap100PTree3)
```

![plot of chunk unnamed-chunk-8](figure/chrs-unnamed-chunk-8-4.png)


#MAKE LENTO PLOTS
#CLADE 1

```r
#get support for each split, and split dot matrix
discoPSuppCl1 <- getSplitSupports(discoPTree1,discoGenos[cl1,])
discoPSplitsCl1 <- getSplitPlot(discoPTree1)

#order by most-supported split
splitSort <- discoPSuppCl1$splits[rev(order(discoPSuppCl1$pro))]
discoPSplitsCl1$splits <- factor(discoPSplitsCl1$splits,levels=splitSort, ordered=T)
discoPSuppCl1$splits <- factor(discoPSuppCl1$splits,levels=splitSort, ordered=T)

label = paste("mean split\nsupport: ",round(mean(discoPSuppCl1[,4]/rowSums(discoPSuppCl1[,4:5])),2),sep="")

discoPSuppP1 <- ggplot(discoPSuppCl1,aes(x=splits)) + geom_bar(aes(y=pro),stat="identity",fill="blue") + geom_bar(aes(y=(anti*-1)),stat="identity",fill="red") + ylim(-1000,1750) + annotate("text",x=5,y=1000,label=label,hjust=0,size=6)
discoSplitP1 <- ggplot(discoPSplitsCl1,aes(x=splits,y=sample,fill=insplit)) + geom_point(shape=21, size=3) + scale_fill_manual(values=c(NA,"black"),guide=F)

grid.arrange(discoPSuppP1, discoSplitP1, ncol=1)
```

![plot of chunk unnamed-chunk-9](figure/chrs-unnamed-chunk-9-1.png)


```r
#get support for each split, and split dot matrix
hap250PSuppCl1 <- getSplitSupports(hap250PTree1,hap250Genos[cl1,])
hap250PSplitsCl1 <- getSplitPlot(hap250PTree1)

#order by most-supported split
splitSort <- hap250PSuppCl1$splits[rev(order(hap250PSuppCl1$pro))]
hap250PSplitsCl1$splits <- factor(hap250PSplitsCl1$splits,levels=splitSort, ordered=T)
hap250PSuppCl1$splits <- factor(hap250PSuppCl1$splits,levels=splitSort, ordered=T)

label = paste("mean split\nsupport: ",round(mean(hap250PSuppCl1[,4]/rowSums(hap250PSuppCl1[,4:5])),2),sep="")


hap250PSuppP1 <- ggplot(hap250PSuppCl1,aes(x=splits)) + geom_bar(aes(y=pro),stat="identity",fill="blue") + geom_bar(aes(y=(anti*-1)),stat="identity",fill="red") + ylim(-1000,1750) + annotate("text",x=5,y=1000,label=label,hjust=0,size=6)
hap250SplitP1 <- ggplot(hap250PSplitsCl1,aes(x=splits,y=sample,fill=insplit)) + geom_point(shape=21, size=3) + scale_fill_manual(values=c(NA,"black"),guide=F)

grid.arrange(hap250PSuppP1, hap250SplitP1, ncol=1)
```

![plot of chunk unnamed-chunk-10](figure/chrs-unnamed-chunk-10-1.png)


```r
#get support for each split, and split dot matrix
hap100PSuppCl1 <- getSplitSupports(hap100PTree1,hap100Genos[cl1,])
hap100PSplitsCl1 <- getSplitPlot(hap100PTree1)

#order by most-supported split
splitSort <- hap100PSuppCl1$splits[rev(order(hap100PSuppCl1$pro))]
hap100PSplitsCl1$splits <- factor(hap100PSplitsCl1$splits,levels=splitSort, ordered=T)
hap100PSuppCl1$splits <- factor(hap100PSuppCl1$splits,levels=splitSort, ordered=T)

label = paste("mean split\nsupport: ",round(mean(hap100PSuppCl1[,4]/rowSums(hap100PSuppCl1[,4:5])),2),sep="")

hap100PSuppP1 <- ggplot(hap100PSuppCl1,aes(x=splits)) + geom_bar(aes(y=pro),stat="identity",fill="blue") + geom_bar(aes(y=(anti*-1)),stat="identity",fill="red") + ylim(-1000,1750) + annotate("text",x=5,y=1000,label=label,hjust=0,size=6)
hap100SplitP1 <- ggplot(hap100PSplitsCl1,aes(x=splits,y=sample,fill=insplit)) + geom_point(shape=21, size=3) + scale_fill_manual(values=c(NA,"black"),guide=F)
grid.arrange(hap100PSuppP1, hap100SplitP1, ncol=1)
```

![plot of chunk unnamed-chunk-11](figure/chrs-unnamed-chunk-11-1.png)


Clade2

```r
#get support for each split, and split dot matrix
discoPSuppCl2 <- getSplitSupports(discoPTree2,discoGenos[cl2,])
discoPSplitsCl2 <- getSplitPlot(discoPTree2)

#order by most-supported split
splitSort <- discoPSuppCl2$splits[rev(order(discoPSuppCl2$pro))]
discoPSplitsCl2$splits <- factor(discoPSplitsCl2$splits,levels=splitSort, ordered=T)
discoPSuppCl2$splits <- factor(discoPSuppCl2$splits,levels=splitSort, ordered=T)

label = paste("mean split\nsupport: ",round(mean(discoPSuppCl2[,4]/rowSums(discoPSuppCl2[,4:5])),2),sep="")

discoSplitP2 <- ggplot(discoPSplitsCl2,aes(x=splits,y=sample,fill=insplit)) + geom_point(shape=21, size=3) + scale_fill_manual(values=c(NA,"black"),guide=F)

discoPSuppP2 <- ggplot(discoPSuppCl2,aes(x=splits)) + geom_bar(aes(y=pro),stat="identity",fill="blue") + geom_bar(aes(y=(anti*-1)),stat="identity",fill="red") + ylim(-1000,1750) + annotate("text",x=2,y=1000,label=label,hjust=0,size=6)

grid.arrange(discoPSuppP2 + blankTheme, discoSplitP2+blankTheme, ncol=1)
```

![plot of chunk unnamed-chunk-12](figure/chrs-unnamed-chunk-12-1.png)


```r
#get support for each split, and split dot matrix
hap250PSuppCl2 <- getSplitSupports(hap250PTree2,hap250Genos[cl2,])
hap250PSplitsCl2 <- getSplitPlot(hap250PTree2)

#order by most-supported split
splitSort <- hap250PSuppCl2$splits[rev(order(hap250PSuppCl2$pro))]
hap250PSplitsCl2$splits <- factor(hap250PSplitsCl2$splits,levels=splitSort, ordered=T)
hap250PSuppCl2$splits <- factor(hap250PSuppCl2$splits,levels=splitSort, ordered=T)

label = paste("mean split\nsupport: ",round(mean(hap250PSuppCl2[,4]/rowSums(hap250PSuppCl2[,4:5])),2),sep="")

hap250SplitP2 <- ggplot(hap250PSplitsCl2,aes(x=splits,y=sample,fill=insplit)) + geom_point(shape=21, size=3) + scale_fill_manual(values=c(NA,"black"),guide=F)

hap250PSuppP2 <- ggplot(hap250PSuppCl2,aes(x=splits)) + geom_bar(aes(y=pro),stat="identity",fill="blue") + geom_bar(aes(y=(anti*-1)),stat="identity",fill="red") + ylim(-1000,1750) + annotate("text",x=2,y=1000,label=label,hjust=0,size=6)
grid.arrange(hap250PSuppP2 + blankTheme, hap250SplitP2+blankTheme, ncol=1)
```

![plot of chunk unnamed-chunk-13](figure/chrs-unnamed-chunk-13-1.png)


```r
#get support for each split, and split dot matrix
hap100PSuppCl2 <- getSplitSupports(hap100PTree2,hap100Genos[cl2,])
hap100PSplitsCl2 <- getSplitPlot(hap100PTree2)

#order by most-supported split
splitSort <- hap100PSuppCl2$splits[rev(order(hap100PSuppCl2$pro))]
hap100PSplitsCl2$splits <- factor(hap100PSplitsCl2$splits,levels=splitSort, ordered=T)
hap100PSuppCl2$splits <- factor(hap100PSuppCl2$splits,levels=splitSort, ordered=T)

label = paste("mean split\nsupport: ",round(mean(hap100PSuppCl2[,4]/rowSums(hap100PSuppCl2[,4:5])),2),sep="")

hap100SplitP2 <- ggplot(hap100PSplitsCl2,aes(x=splits,y=sample,fill=insplit)) + geom_point(shape=21, size=3) + scale_fill_manual(values=c(NA,"black"),guide=F)
hap100PSuppP2 <- ggplot(hap100PSuppCl2,aes(x=splits)) + geom_bar(aes(y=pro),stat="identity",fill="blue") + geom_bar(aes(y=(anti*-1)),stat="identity",fill="red") + ylim(-1000,1750) + annotate("text",x=2,y=1000,label=label,hjust=0,size=6)
grid.arrange(hap100PSuppP2 + blankTheme, hap100SplitP2+blankTheme, ncol=1)
```

![plot of chunk unnamed-chunk-14](figure/chrs-unnamed-chunk-14-1.png)


Clade3

```r
#get support for each split, and split dot matrix
discoPSuppCl3 <- getSplitSupports(discoPTree3,discoGenos[cl3,])
discoPSplitsCl3 <- getSplitPlot(discoPTree3)

#order by most-supported split
splitSort <- discoPSuppCl3$splits[rev(order(discoPSuppCl3$pro))]
discoPSplitsCl3$splits <- factor(discoPSplitsCl3$splits,levels=splitSort, ordered=T)
discoPSuppCl3$splits <- factor(discoPSuppCl3$splits,levels=splitSort, ordered=T)

discoSplitP3 <- ggplot(discoPSplitsCl3,aes(x=splits,y=sample,fill=insplit)) + geom_point(shape=21, size=3) + scale_fill_manual(values=c(NA,"black"),guide=F)

discoPSuppP3 <- ggplot(discoPSuppCl3,aes(x=splits)) + geom_bar(aes(y=pro),stat="identity",fill="blue") + geom_bar(aes(y=(anti*-1)),stat="identity",fill="red") + ylim(-1000,1750)

grid.arrange(discoPSuppP3 + blankTheme, discoSplitP3+blankTheme, ncol=1)
```

![plot of chunk unnamed-chunk-15](figure/chrs-unnamed-chunk-15-1.png)


```r
#get support for each split, and split dot matrix
hap250PSuppCl3 <- getSplitSupports(hap250PTree3,hap250Genos[cl3,])
hap250PSplitsCl3 <- getSplitPlot(hap250PTree3)

#order by most-supported split
splitSort <- hap250PSuppCl3$splits[rev(order(hap250PSuppCl3$pro))]
hap250PSplitsCl3$splits <- factor(hap250PSplitsCl3$splits,levels=splitSort, ordered=T)
hap250PSuppCl3$splits <- factor(hap250PSuppCl3$splits,levels=splitSort, ordered=T)

hap250SplitP3 <- ggplot(hap250PSplitsCl3,aes(x=splits,y=sample,fill=insplit)) + geom_point(shape=21, size=3) + scale_fill_manual(values=c(NA,"black"),guide=F)

hap250PSuppP3 <- ggplot(hap250PSuppCl3,aes(x=splits)) + geom_bar(aes(y=pro),stat="identity",fill="blue") + geom_bar(aes(y=(anti*-1)),stat="identity",fill="red") + ylim(-1000,1750)
grid.arrange(hap250PSuppP3 + blankTheme, hap250SplitP3+blankTheme, ncol=1)
```

![plot of chunk unnamed-chunk-16](figure/chrs-unnamed-chunk-16-1.png)


```r
#get support for each split, and split dot matrix
hap100PSuppCl3 <- getSplitSupports(hap100PTree3,hap100Genos[cl3,])
hap100PSplitsCl3 <- getSplitPlot(hap100PTree3)

#order by most-supported split
splitSort <- hap100PSuppCl3$splits[rev(order(hap100PSuppCl3$pro))]
hap100PSplitsCl3$splits <- factor(hap100PSplitsCl3$splits,levels=splitSort, ordered=T)
hap100PSuppCl3$splits <- factor(hap100PSuppCl3$splits,levels=splitSort, ordered=T)

hap100SplitP3 <- ggplot(hap100PSplitsCl3,aes(x=splits,y=sample,fill=insplit)) + geom_point(shape=21, size=3) + scale_fill_manual(values=c(NA,"black"),guide=F)
hap100PSuppP3 <- ggplot(hap100PSuppCl3,aes(x=splits)) + geom_bar(aes(y=pro),stat="identity",fill="blue") + geom_bar(aes(y=(anti*-1)),stat="identity",fill="red") + ylim(-1000,1750)
grid.arrange(hap100PSuppP3 + blankTheme, hap100SplitP3+blankTheme, ncol=1)
```

![plot of chunk unnamed-chunk-17](figure/chrs-unnamed-chunk-17-1.png)






```r
opts_chunk$set(fig.width=21, fig.height=12)
opts_chunk$set(dev=c('png','postscript'))
#opts_chunk$set(dev=c('png'))
```


```r
blank <- grid.rect(gp=gpar(col="white"))
```

![plot of chunk unnamed-chunk-19](figure/chrs-unnamed-chunk-19-1.png)

```r
grid.arrange(blank,discoPSuppP1 + blankTheme, 
             blank,hap250PSuppP1 + blankTheme,
             blank,hap100PSuppP1 + blankTheme, 
             
             blank,discoSplitP1+blankTheme, 
             blank,hap250SplitP1+blankTheme, 
             blank,hap100SplitP1+blankTheme, 
             ncol=6,widths=c(2,6,2,6,2,6))
```

![plot of chunk unnamed-chunk-19](figure/chrs-unnamed-chunk-19-2.png)



```r
opts_chunk$set(fig.width=12, fig.height=8)
opts_chunk$set(dev=c('png','postscript'))
#opts_chunk$set(dev=c('png'))
```


```r
grid.arrange(blank,discoPSuppP2 + blankTheme, 
             blank,hap250PSuppP2 + blankTheme,
             blank,hap100PSuppP2 + blankTheme, 
             
             blank,discoSplitP2+blankTheme, 
             blank,hap250SplitP2+blankTheme, 
             blank,hap100SplitP2+blankTheme, 
             ncol=6,widths=c(2,6,2,6,2,6))
```

![plot of chunk unnamed-chunk-21](figure/chrs-unnamed-chunk-21-1.png)



```r
opts_chunk$set(fig.width=8, fig.height=8)
opts_chunk$set(dev=c('png','postscript'))
#opts_chunk$set(dev=c('png'))
```


```r
grid.arrange(blank,discoPSuppP3 + blankTheme, 
             blank,hap250PSuppP3 + blankTheme,
             blank,hap100PSuppP3 + blankTheme, 
             
             blank,discoSplitP3+blankTheme, 
             blank,hap250SplitP3+blankTheme, 
             blank,hap100SplitP3+blankTheme, 
             ncol=6,widths=c(2,6,2,6,2,6))
```

![plot of chunk unnamed-chunk-23](figure/chrs-unnamed-chunk-23-1.png)



```r
max(c(c(discoPSuppCl1$anti,discoPSuppCl2$anti,discoPSuppCl3$anti),
    c(hap100PSuppCl1$anti,hap100PSuppCl2$anti,hap100PSuppCl3$anti),
    c(hap250PSuppCl1$anti,hap250PSuppCl2$anti,hap250PSuppCl3$anti)))
```

```
## [1] 769
```

```r
max(c(c(discoPSuppCl1$pro,discoPSuppCl2$pro,discoPSuppCl3$pro),
    c(hap100PSuppCl1$pro,hap100PSuppCl2$pro,hap100PSuppCl3$pro),
    c(hap250PSuppCl1$pro,hap250PSuppCl2$pro,hap250PSuppCl3$pro)))
```

```
## [1] 1559
```

```r
#percent support Discovar:
write.table(round(discoPSuppCl1[,3:5] / rowSums(discoPSuppCl1[,3:5]),2),sep="\t",quote=F)
```

```
## irrelevant	pro	anti
## 1	0.69	0.31	0
## 2	0.7	0.26	0.03
## 3	0.92	0.05	0.04
## 4	0.9	0.07	0.03
## 5	0.91	0.06	0.03
## 6	0.95	0.04	0.02
## 7	0.9	0.07	0.02
## 8	0.94	0.05	0.01
```

```r
mean(discoPSuppCl1[,4]/rowSums(discoPSuppCl1[,4:5]))
```

```
## [1] 0.7529791
```

```r
write.table(discoPSuppCl1[,3:5],quote=F) 
```

```
## irrelevant pro anti
## 1 1053 471 0
## 2 1155 433 57
## 3 1512 75 59
## 4 1482 109 50
## 5 1453 98 50
## 6 1463 58 27
## 7 1460 116 40
## 8 1475 71 17
```

```r
#percent support Hap250:
write.table(round(hap250PSuppCl1[,3:5] / rowSums(hap250PSuppCl1[,3:5]),2),sep="\t",quote=F)
```

```
## irrelevant	pro	anti
## 1	0.54	0.23	0.23
## 2	0.75	0.01	0.24
## 3	0.76	0.01	0.23
## 4	0.81	0.01	0.18
## 5	0.87	0.03	0.1
## 6	0.54	0.23	0.23
## 7	0.79	0.04	0.17
## 8	0.86	0.03	0.12
```

```r
mean(hap250PSuppCl1[,4]/rowSums(hap250PSuppCl1[,4:5]))
```

```
## [1] 0.2099631
```

```r
write.table(hap250PSuppCl1[,3:5],quote=F)
```

```
## irrelevant pro anti
## 1 1415 600 602
## 2 1959 19 639
## 3 1989 18 610
## 4 2123 23 471
## 5 2282 66 267
## 6 1415 600 602
## 7 2064 112 437
## 8 2236 66 311
```

```r
#percent support Hap100:
write.table(round(hap100PSuppCl1[,3:5] / rowSums(hap100PSuppCl1[,3:5]),2),sep="\t",quote=F)
```

```
## irrelevant	pro	anti
## 1	0.53	0.47	0
## 2	0.65	0.13	0.22
## 3	0.77	0.01	0.22
## 4	0.8	0.01	0.2
## 5	0.88	0.01	0.11
## 6	0.89	0	0.11
## 7	0.82	0.03	0.15
## 8	0.89	0.01	0.11
```

```r
mean(hap100PSuppCl1[,4]/rowSums(hap100PSuppCl1[,4:5]))
```

```
## [1] 0.2208524
```

```r
write.table(hap100PSuppCl1[,3:5],quote=F)
```

```
## irrelevant pro anti
## 1 1733 1559 0
## 2 2294 470 769
## 3 2734 33 766
## 4 2813 19 696
## 5 3112 21 394
## 6 3128 16 388
## 7 2908 109 516
## 8 3132 22 378
```



```r
#clade 2/3 pc disco support
write.table(round(discoPSuppCl2[,3:5] / rowSums(discoPSuppCl2[,3:5]),2),sep="\t",quote=F)
```

```
## irrelevant	pro	anti
## 1	0.68	0.32	0
## 2	0.8	0.16	0.04
## 3	0.85	0.12	0.03
```

```r
write.table(round(discoPSuppCl3[,3:5] / rowSums(discoPSuppCl3[,3:5]),2),sep="\t",quote=F)
```

```
## irrelevant	pro	anti
## 1	0.5	0.5	0
```

```r
#clade 2/3 pc HC250 support
write.table(round(hap250PSuppCl2[,3:5] / rowSums(hap250PSuppCl2[,3:5]),2),sep="\t",quote=F)
```

```
## irrelevant	pro	anti
## 1	0.83	0.17	0
## 2	0.69	0.15	0.15
## 3	0.7	0.1	0.2
```

```r
write.table(round(hap250PSuppCl3[,3:5] / rowSums(hap250PSuppCl3[,3:5]),2),sep="\t",quote=F)
```

```
## irrelevant	pro	anti
## 1	0.49	0.51	0
```

```r
#clade 2/3 pc HC100 support
write.table(round(hap100PSuppCl2[,3:5] / rowSums(hap100PSuppCl2[,3:5]),2),sep="\t",quote=F)
```

```
## irrelevant	pro	anti
## 1	0.81	0.19	0
## 2	0.67	0.11	0.22
## 3	0.7	0.16	0.14
```

```r
write.table(round(hap100PSuppCl3[,3:5] / rowSums(hap100PSuppCl3[,3:5]),2),sep="\t",quote=F)
```

```
## irrelevant	pro	anti
## 1	0.59	0.41	0
```


library(ggplot2)
library(reshape2)
library(knitr)

```r
opts_chunk$set(fig.width=13, fig.height=12)
```


```r
getGenoTab <- function(filename) {
  genotab <- read.table(filename,sep="\t",header=T,stringsAsFactors=F,na.strings=c('.','. '))
  genotab <- unique(genotab)
  genotab$I <- 1:length(genotab$chr)
  #remove ref calls
  #genotab <- genotab[rowSums(genotab[,5:dim(genotab)[2]],na.rm=T)!=0,]
  
  #colnames(genotab)[1:2]<-c("chr","pos")
  genotab.m <- melt(genotab,id.vars =c("I","chr","pos","type","subtype","alleles"),variable.name ="indv")
  genotab.m$value <- as.factor(genotab.m$value)
  #sort levels
  genotab.m$indv <- factor(genotab.m$indv,levels=sort(levels(genotab.m$indv)))
  
  genotab.m <- genotab.m[!is.na(genotab.m$value),]
  genotab.m$value <- as.numeric(as.character(genotab.m$value))
  
  hNo=0
  vlast=0
  genotab.m$haplo=0
  for (i in c(1:dim(genotab.m)[[1]])) {
    if (genotab.m[i,"value"]==0) {haplo=0}
    else if (genotab.m[i,"value"]>=1) {
      if (vlast==0) {hNo <- hNo+1}
        genotab.m[i,"haplo"]<-hNo  
      }
    vlast=genotab.m[i,"value"]
  }
  genotab.m$value <- factor(genotab.m$value)
  genotab.m$haplo <- factor(genotab.m$haplo)
  genotab.m
}
```


```r
networkorder <- c("Th134.11","Th117.11","Th106.11","Th086.07","Th106.09","Th196.12","Th074.13","Th132.11","Th162.12","Th230.12")

clade3F <- "Thies_all_manual.PASS.Cls.miss0.5.LMRG.HAP.DISCORDS.clade3.alleles.tab.txt"
genotab3 <- getGenoTab(clade3F)
genotab3$indv <- factor(genotab3$indv,levels = networkorder,ordered=T)

clade3NS <- "Thies_all_manual.PASS.Cls.miss0.5.LMRG.HAP.DISCORDS.clade3.NOSNGL.alleles.tab.txt"
genotab3NS <- getGenoTab(clade3NS)
genotab3NS$indv <- factor(genotab3NS$indv,levels = networkorder,ordered=T)
```


```r
cols=c("grey","light green","dark green","orange","red","dark red")
names(cols)<-c(0,1,2,3,4,5)
colsc <- scale_color_manual(values = cols)
```




```r
#Clade3

#all vars
ggplot(genotab3,aes(x=pos,y=indv,colour=value,group=haplo,shape=type)) + geom_point(size=3) + 
  geom_line(data=subset(genotab3,value!=0)) + facet_grid(chr ~ .) + xlim(0,3.5e6) +colsc
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-1.png)

```r
#SNPs only
ggplot(subset(genotab3,type=="SNP"),aes(x=pos,y=indv,colour=value,group=haplo,shape=type)) + geom_point(alpha=0.6,size=3) + 
  geom_line(data=subset(genotab3,value!=0 & type=="SNP")) + facet_grid(chr ~ type) + xlim(0,3.5e6) +colsc
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-2.png)

```r
#INDELs only
ggplot(subset(genotab3,type=="INDEL"),aes(x=pos,y=indv,colour=value,group=haplo,shape=type)) + geom_point(alpha=0.6,size=3) + 
  geom_line(data=subset(genotab3,value!=0 & type=="INDEL")) + facet_grid(chr ~ type) + xlim(0,3.5e6) +colsc
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-3.png)




```r
#Clade3 NO SINGLETONS

#all vars
ggplot(genotab3NS,aes(x=pos,y=indv,colour=value,group=haplo,shape=type)) + geom_point(size=3) + 
  geom_line(data=subset(genotab3NS,value!=0)) + facet_grid(chr ~ .) + xlim(0,3.5e6) +colsc
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-1.png)

```r
#SNPs only
ggplot(subset(genotab3NS,type=="SNP"),aes(x=pos,y=indv,colour=value,group=haplo,shape=type)) + geom_point(alpha=0.6,size=3) + 
  geom_line(data=subset(genotab3NS,value!=0 & type=="SNP")) + facet_grid(chr ~ type) + xlim(0,3.5e6) +colsc
```

```
## geom_path: Each group consists of only one observation. Do you need to
## adjust the group aesthetic?
## geom_path: Each group consists of only one observation. Do you need to
## adjust the group aesthetic?
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-2.png)

```r
#INDELs only
ggplot(subset(genotab3NS,type=="INDEL"),aes(x=pos,y=indv,colour=value,group=haplo,shape=type)) + geom_point(alpha=0.6,size=3) + 
  geom_line(data=subset(genotab3NS,value!=0 & type=="INDEL")) + facet_grid(chr ~ type) + xlim(0,3.5e6) +colsc
```

```
## geom_path: Each group consists of only one observation. Do you need to
## adjust the group aesthetic?
## geom_path: Each group consists of only one observation. Do you need to
## adjust the group aesthetic?
## geom_path: Each group consists of only one observation. Do you need to
## adjust the group aesthetic?
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-3.png)


```r
#Clade3 NO SINGLETONS INDEXED

#SNPs only
ggplot(genotab3NS,aes(x=I,y=indv,colour=value,group=pos,shape=type)) + geom_point(alpha=0.6,size=3) + 
  geom_line(data=subset(genotab3NS,value != 0)) + facet_grid(type ~ .)  +colsc
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6-1.png)

```r
#INDELs only
```
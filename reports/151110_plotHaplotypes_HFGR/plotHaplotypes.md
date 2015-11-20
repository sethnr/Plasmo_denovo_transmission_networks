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
  
  #remove ref calls
  #genotab <- genotab[rowSums(genotab[,5:dim(genotab)[2]],na.rm=T)!=0,]
  
  #colnames(genotab)[1:2]<-c("chr","pos")
  genotab.m <- melt(genotab,id.vars =c("chr","pos","type","alleles"),variable.name ="indv")
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
clade1F <- "HFGR_all.PASS.HFGR.miss0.5.LMRG.HAP.DISCORD.alleles.tab.txt"
genotab1 <- getGenoTab(clade1F)

newOrder <- c( "HFGRII.2.10x","HFGRII.1.30x"  , "HFGRII.1.200x" , "HFGRII.2.200x" , "HFGRII.2.300x"  ,
"HFGRIII.1.30x",  "HFGRIII.1.60x",  "HFGRIII.2.60x",  "HFGRIII.3.200x")
newOrder <- rev(newOrder)
genotab1$indv <- factor(genotab1$indv,levels = newOrder,ordered=T)
```


```r
cols=c("grey","dark green","orange","dark red")
names(cols)<-c(0,1,2,3)
colsc <- scale_color_manual(values = cols)
```



```r
#HFGR - all

#all vars
ggplot(genotab1,aes(x=pos,y=indv,colour=value,group=haplo,shape=type)) + geom_point(size=3) + 
  geom_line(data=subset(genotab1,value!=0)) + facet_grid(chr ~ .) + xlim(0,3.5e6) +colsc
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-1.png) 

```r
#SNPs only
ggplot(subset(genotab1,type=="SNP"),aes(x=pos,y=indv,colour=value,group=haplo,shape=type)) + geom_point(alpha=0.6,size=3) + 
  geom_line(data=subset(genotab1,value!=0 & type=="SNP")) + facet_grid(chr ~ type) + xlim(0,3.5e6) +colsc
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-2.png) 

```r
#INDELs only
ggplot(subset(genotab1,type=="INDEL"),aes(x=pos,y=indv,colour=value,group=haplo,shape=type)) + geom_point(alpha=0.6,size=3) + 
  geom_line(data=subset(genotab1,value!=0 & type=="INDEL")) + facet_grid(chr ~ type) + xlim(0,3.5e6) +colsc
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-3.png) 

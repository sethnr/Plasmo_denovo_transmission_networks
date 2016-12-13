library(ggplot2)
library(reshape2)
library(knitr)
library(igraph)


```r
opts_chunk$set(fig.width=13, fig.height=12)
```


```r
getGenoTab <- function(filename) {
  genotab <- read.table(filename,sep="\t",header=T,stringsAsFactors=F,na.strings=c('.','. '))
#  genotab <- unique(genotab)
  genotab$I <- 1:length(genotab$chr)
  #remove ref calls
  #genotab <- genotab[rowSums(genotab[,5:dim(genotab)[2]],na.rm=T)!=0,]
  
  #colnames(genotab)[1:2]<-c("chr","pos")
  genotab.m <- melt(genotab,id.vars =c("I","chr","pos","type","subtype","alleles"),variable.name ="indv")
  #genotab.m$value <- as.factor(genotab.m$value)
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
  #genotab.m$value <- factor(genotab.m$value)
  genotab.m$haplo <- factor(genotab.m$haplo)
  genotab.m
}
```


```r
#networkorder <- c("Th134.11","Th117.11","Th106.11","Th086.07","Th106.09","Th196.12","Th074.13","Th132.11","Th162.12","Th230.12")
INDELnetworkorder <- c("Th134.11","Th117.11","Th106.11","Th086.07","Th106.09","Th196.12","Th074.13","Th132.11","Th162.12","Th230.12")
SNPnetworkorder <- c("Th086.07","Th106.09","Th106.11","Th117.11","Th134.11","Th230.12","Th162.12","Th074.13","Th132.11","Th196.12")

clade3F <- "clade29.alleles.tab.txt"
genotab3 <- getGenoTab(clade3F)
genotab3$indv <- factor(genotab3$indv,levels = networkorder,ordered=T)
```

```
## Error in factor(genotab3$indv, levels = networkorder, ordered = T): object 'networkorder' not found
```

```r
colnames(genotab3)[[8]]<-"allele"

clade3NSF <- "clade29.alleles.NOSNGL.tab.txt"
genotab3NS <- getGenoTab(clade3NSF)
genotab3NS$indv <- factor(genotab3NS$indv,levels = networkorder,ordered=T)
```

```
## Error in factor(genotab3NS$indv, levels = networkorder, ordered = T): object 'networkorder' not found
```

```r
colnames(genotab3NS)[[8]]<-"allele"


clade3HPSF <- "clade29.SNP.graphml.HPLSY.txt"
genotab3HPS <- getGenoTab(clade3HPSF)
genotab3HPS$indv <- factor(genotab3HPS$indv,levels = networkorder,ordered=T)
```

```
## Error in factor(genotab3HPS$indv, levels = networkorder, ordered = T): object 'networkorder' not found
```

```r
colnames(genotab3HPS)[[8]]<-"HP_SNP"

clade3HPIF <- "clade29.INDEL.graphml.HPLSY.txt"
genotab3HPI <- getGenoTab(clade3HPIF)
genotab3HPI$indv <- factor(genotab3HPI$indv,levels = networkorder,ordered=T)
```

```
## Error in factor(genotab3HPI$indv, levels = networkorder, ordered = T): object 'networkorder' not found
```

```r
colnames(genotab3HPI)[[8]]<-"HP_INDEL"

genotab3HPS <- genotab3HPS[,c(2:3,7:8)]
genotab3HPI <- genotab3HPI[,c(2:3,7:8)]
dim(merge(genotab3HPS,genotab3HPI,by=c("chr","pos","indv"),all=T))
```

```
## [1] 13690     5
```

```r
#dim(merge(genotab3,merge(genotab3HPS,genotab3HPI,by=c("chr","pos","indv"),all=T)))
genotab <- merge(genotab3,merge(genotab3HPS,genotab3HPI,by=c("chr","pos","indv"),all=T))
genotabNS <- merge(genotab3NS,merge(genotab3HPS,genotab3HPI,by=c("chr","pos","indv"),all=T))
dim(genotab)
```

```
## [1] 13392    11
```

```r
dim(genotabNS)
```

```
## [1] 4777   11
```

```r
dim(genotab3HPS)
```

```
## [1] 13570     4
```

```r
dim(genotab3HPI)
```

```
## [1] 13570     4
```

```r
genotabNS$indv <- factor(genotabNS$indv,levels = networkorder,ordered=T)
```

```
## Error in factor(genotabNS$indv, levels = networkorder, ordered = T): object 'networkorder' not found
```

```r
genotab$indv <- factor(genotab$indv,levels = networkorder,ordered=T)
```

```
## Error in factor(genotab$indv, levels = networkorder, ordered = T): object 'networkorder' not found
```

```r
genotab$allele <- factor(genotab$allele)
genotabNS$allele <- factor(genotabNS$allele)
```


```r
cols=c("grey","light green","dark green","orange","red","dark red")
names(cols)<-c(0,1,2,3,4,5)
colsc <- scale_color_manual(values = cols)

fills=c("#BBBBBB","#DDDDDD")
names(fills)<-c("0","1")
fillsc <- scale_fill_manual(values = fills)
```




```r
#Clade3

#all vars
ggplot(genotab,aes(x=pos,y=indv,colour=allele,group=haplo,shape=type)) + geom_point(size=3) + 
  geom_line(data=subset(genotab,allele!=0)) + facet_grid(chr ~ .) + xlim(0,3.5e6) +colsc
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-1.png)

```r
#SNP net homoplasies only
#ggplot(genotab,aes(x=pos,y=indv,colour=HP_SNP,group=haplo,shape=type)) + geom_point(alpha=0.6,size=3) + 
#  geom_line(data=subset(genotab,allele!=0)) + facet_grid(chr ~ type) + xlim(0,3.5e6) +colsc
ggplot(genotab,aes(x=pos,y=indv,colour=allele,group=HP_SNP,shape=type)) + geom_point(alpha=0.6,size=3) + 
  geom_line(data=subset(genotab,HP_SNP!=0),colour="red") + facet_grid(chr ~ type) + xlim(0,3.5e6) +colsc
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-2.png)

```r
#INDEL net homoplasies only
ggplot(genotab,aes(x=pos,y=indv,colour=allele,group=HP_INDEL,shape=type)) + geom_point(alpha=0.6,size=3) + 
  geom_line(data=subset(genotab,HP_INDEL!=0),color="red") + facet_grid(chr ~ type) + xlim(0,3.5e6) +colsc
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-3.png)


```r
opts_chunk$set(dev=c('png','pdf')) 
```


```r
#Clade3 no singletons

#all vars
ggplot(genotabNS,aes(x=pos,y=indv,colour=allele,group=haplo,shape=type)) + geom_point(size=3) + 
  geom_line(data=subset(genotabNS,allele!=0)) + facet_grid(chr ~ .) + xlim(0,3.5e6) +colsc
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6-1.png)

```r
#SNP net homoplasies only
ggplot(genotabNS,aes(x=pos,y=indv,colour=allele,group=HP_SNP,shape=type)) + geom_point(alpha=0.6,size=3) + 
  geom_line(data=subset(genotabNS,HP_SNP!=0),colour="red") + facet_grid(chr ~ type) + xlim(0,3.5e6) +colsc
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6-2.png)

```r
#INDEL net homoplasies only
ggplot(genotabNS,aes(x=pos,y=indv,colour=allele,group=HP_INDEL,shape=type)) + geom_point(alpha=0.6,size=3) + 
  geom_line(data=subset(genotabNS,HP_INDEL!=0),colour="red") + facet_grid(chr ~ type) + xlim(0,3.5e6) +colsc
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6-3.png)



```r
#Clade3 NO SINGLETONS INDEXED

#SNP net homoplasies
genotabNS$indv <- factor(genotabNS$indv,levels = SNPnetworkorder,ordered=T)
ggplot(genotabNS,aes(x=I,y=indv,colour=allele,group=HP_SNP,shape=type)) + 
  geom_point(alpha=0.6,size=3) +geom_line(data=subset(genotabNS,HP_SNP>0),colour="red") + facet_grid(type ~ .)  +colsc
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-1.png)

```r
#INDEL net homoplasies
genotabNS$indv <- factor(genotabNS$indv,levels = INDELnetworkorder,ordered=T)
ggplot(genotabNS,aes(x=I,y=indv,colour=allele,group=HP_INDEL,shape=type)) + 
  geom_point(alpha=0.6,size=3) + geom_line(data=subset(genotabNS,HP_INDEL >0),colour="red") + facet_grid(type ~ .)  +colsc
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-2.png)


```r
chroms <- cbind(merge(aggregate(genotab$I,by = list(genotab$chr),FUN=min),aggregate(genotab$I,by = list(genotab$chr),FUN=max),by="Group.1"),rep(c("0","1"),7))
colnames(chroms) <- c("chr","min","max","col")
chromsNS <- cbind(merge(aggregate(genotabNS$I,by = list(genotabNS$chr),FUN=min),aggregate(genotabNS$I,by = list(genotabNS$chr),FUN=max),by="Group.1"),rep(c("0","1"),7))
colnames(chromsNS) <- c("chr","min","max","col")


#both var types together
genotabNS$indv <- factor(genotabNS$indv,levels = SNPnetworkorder,ordered=T)
ggplot(genotabNS,aes(x=I,y=indv,group=HP_SNP,shape=type)) + 
  geom_point(aes(colour=allele),alpha=0.6,size=3) +colsc +geom_line(data=subset(genotabNS,HP_SNP>0),colour="red") + 
  geom_rect(data=chromsNS,aes(xmin=min,xmax=max,ymin=0,ymax=0.5,fill=col),inherit.aes=F,colour="#BBBBBB") + fillsc
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-1.png)

```r
genotabNS$indv <- factor(genotabNS$indv,levels = INDELnetworkorder,ordered=T)
ggplot(genotabNS,aes(x=I,y=indv,colour=allele,group=HP_INDEL,shape=type)) + 
  geom_point(alpha=0.6,size=3) +colsc +geom_line(data=subset(genotabNS,HP_INDEL >0),colour="red") + 
  geom_rect(data=chromsNS,aes(xmin=min,xmax=max,ymin=0,ymax=0.5,fill=col),inherit.aes=F,colour="#BBBBBB") + fillsc
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-2.png)

```r
#inc singletons
genotab$indv <- factor(genotab$indv,levels = SNPnetworkorder,ordered=T)
ggplot(genotab,aes(x=I,y=indv,group=HP_SNP,shape=type)) + 
  geom_point(aes(colour=allele),alpha=0.6,size=3) +colsc +geom_line(data=subset(genotab,HP_SNP>0),colour="red") + 
  geom_rect(data=chroms,aes(xmin=min,xmax=max,ymin=0,ymax=0.5,fill=col),inherit.aes=F,colour="#BBBBBB") + fillsc
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-3.png)

```r
genotab$indv <- factor(genotab$indv,levels = INDELnetworkorder,ordered=T)
ggplot(genotab,aes(x=I,y=indv,group=HP_INDEL,shape=type)) + 
  geom_point(aes(colour=allele),alpha=0.6,size=3) +colsc + geom_line(data=subset(genotab,HP_INDEL>0),colour="red") + 
  geom_rect(data=chroms,aes(xmin=min,xmax=max,ymin=0,ymax=0.5,fill=col),inherit.aes=F,colour="#BBBBBB") + fillsc
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-4.png)



```r
opts_chunk$set(fig.width=9, fig.height=9,dev=c('png','postscript'))
```


```r
#SNP net homoplasies
snpnet <- read.graph("clade29.SNP.graphml.HPLSY.xml",format="graphml")
sfrlay <- layout.fruchterman.reingold(snpnet)
plot(snpnet,layout=sfrlay,edge.color="black",edge.label=E(snpnet)$homoplasy_count,edge.label.cex=1.5,edge.label.family="Arial",vertex.label.family="Arial",vertex.size=23)
```

```
## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database
```

```
## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10-1.png)

```r
plot(snpnet,layout=sfrlay,edge.color="black",edge.label=round(E(snpnet)$homoplasy_percent,2),edge.label.cex=1.5,edge.label.family="Arial",vertex.label.family="Arial",vertex.size=23)
```

```
## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database
```

```
## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database
```

```
## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10-2.png)

```r
#INDEL net homoplasies
indelnet <- read.graph("clade29.INDEL.graphml.HPLSY.xml",format="graphml")
ifrlay <- layout.fruchterman.reingold(indelnet)
plot(indelnet,layout=ifrlay,edge.color="black",edge.label=E(indelnet)$homoplasy_count,edge.label.cex=1.5,edge.label.family="Arial",vertex.label.family="Arial",vertex.size=23)
```

```
## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database
```

```
## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database
```

```
## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10-3.png)

```r
plot(indelnet,layout=ifrlay,edge.color="black",edge.label=round(E(indelnet)$homoplasy_percent,2),edge.label.cex=1.5,edge.label.family="Arial",vertex.label.family="Arial",vertex.size=23)
```

```
## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database
```

```
## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database
```

```
## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10-4.png)



```r
#SNP net homoplasies
snpnet <- read.graph("clade29.SNP.graphml.HPLSY.xml",format="graphml")
streelay <- layout_as_tree(snpnet,flip.y = F)[,c(2,1)]
plot(snpnet,layout=streelay,edge.color="black",edge.label=E(snpnet)$homoplasy_count,edge.label.cex=1.5,edge.label.family="Arial",vertex.label.family="Arial",vertex.size=23)
```

```
## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database
```

```
## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-11-1.png)

```r
plot(snpnet,layout=streelay,edge.color="black",edge.label=round(E(snpnet)$homoplasy_percent,2),edge.label.cex=1.5,edge.label.family="Arial",vertex.label.family="Arial",vertex.size=23)
```

```
## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database
```

```
## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database
```

```
## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-11-2.png)

```r
#INDEL net homoplasies
indelnet <- read.graph("clade29.INDEL.graphml.HPLSY.xml",format="graphml")
itreelay <- layout_as_tree(indelnet,flip.y = F)[,c(2,1)]
plot(indelnet,layout=itreelay,edge.color="black",edge.label=E(indelnet)$homoplasy_count,edge.label.cex=1.5,edge.label.family="Arial",vertex.label.family="Arial",vertex.size=23)
```

```
## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database
```

```
## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database
```

```
## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-11-3.png)

```r
plot(indelnet,layout=itreelay,edge.color="black",edge.label=round(E(indelnet)$homoplasy_percent,2),edge.label.cex=1.5,edge.label.family="Arial",vertex.label.family="Arial",vertex.size=23)
```

```
## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database
```

```
## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database
```

```
## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-11-4.png)


```r
#lay <- layout_as_tree(snpnet,flip.y = T)
slay <- cbind(as.numeric(as.factor(V(indelnet)$year)),
              c(0,1,2,3,8,4,6,9,5,7))
plot(snpnet,layout=slay,edge.color="black",edge.label.family="Arial",vertex.label.family="Arial",vertex.shape="rectangle",vertex.size=25,vertex.size2=8)
```

```
## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database
```

```
## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12-1.png)

```r
ilay <- cbind(as.numeric(as.factor(V(indelnet)$year)),
              c(0,1,-1,-2,4,-3,5,2,6,3))
plot(indelnet,layout=ilay,edge.color="black",edge.label.family="Arial",vertex.label.family="Arial",vertex.shape="rectangle",vertex.size=25,vertex.size2=8)
```

```
## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database
```

```
## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database

## Warning in text.default(lc.x, lc.y, labels = edge.labels, col =
## edge.label.color, : font family 'Arial' not found in PostScript font
## database
```

```
## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database

## Warning in text.default(x, y, labels = labels, col = label.color, family =
## label.family, : font family 'Arial' not found in PostScript font database
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12-2.png)
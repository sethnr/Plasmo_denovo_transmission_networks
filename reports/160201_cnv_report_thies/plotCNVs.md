library(ggplot2)
library(reshape2)
library(knitr)

```r
opts_chunk$set(fig.width=13, fig.height=12)

header <- c("CHR","ST","EN",	
"SM-7LV7V", "SM-7LV7W",	"SM-7LV7X",	"SM-7LV7Y",	
"SM-7LV7Z",	"SM-7LV81",	"SM-7LV82",	"SM-7LV83",	
"SM-7LV84",	"SM-7LV85",	"SM-7LV86",	"SM-7LV87",	
"SM-7LV88",	"SM-7LV89",	"SM-7LV8A",	"SM-7LV8B",	
"SM-7LV8C",	"SM-7LV8D")


header <- c("CHR","ST","EN",
            "Th086.07","Th106.09","Th106.11","Th117.11",
            "Th132.11","Th134.11","Th162.12","Th196.12",
            "Th230.12","Th074.13","Th166.12","Th092.13",
            "Th211.13","Th245.13","Th246.13","Th068.12",
            "Th061.13","Th095.13")

tree_order <- rev(c("Th068.12","Th061.13","Th095.13",
"Th086.07","Th117.11","Th134.11","Th106.11","Th162.12","Th230.12","Th074.13","Th106.09","Th132.11","Th196.12",
"Th245.13","Th246.13","Th211.13","Th166.12","Th092.13"))

c1<-c("Th068.12","Th061.13","Th095.13")
c2<-c("Th086.07","Th117.11","Th134.11","Th106.11","Th162.12","Th230.12","Th074.13","Th106.09","Th132.11","Th196.12")
c3<-c("Th245.13","Th246.13","Th211.13","Th166.12","Th092.13")
```


```r
getGenoTab <- function(filename) {
  genotab <- read.table(filename,sep="\t",header=T,stringsAsFactors=F,na.strings=c('.','. '))
  genotab <- unique(genotab)
  
  #remove ref calls
  colnames(genotab) <- header
  genotab.m <- melt(genotab,id.vars =c("CHR","ST","EN"),variable.name ="indv")
  genotab.m$value <- as.factor(genotab.m$value)
  #sort levels
  genotab.m$indv <- factor(genotab.m$indv,levels=sort(levels(genotab.m$indv)))
  
  genotab.m <- genotab.m[!is.na(genotab.m$value),]
  
#   hNo=0
#   vlast=0
#   genotab.m$haplo=0
#   for (i in c(1:dim(genotab.m)[[1]])) {
#     if (genotab.m[i,"value"]==0) {haplo=0}
#     else if (genotab.m[i,"value"]>=1) {
#       if (vlast==0) {hNo <- hNo+1}
#         genotab.m[i,"haplo"]<-hNo  
#       }
#     vlast=genotab.m[i,"value"]
#   }
  genotab.m$value <- factor(genotab.m$value)
#   genotab.m$haplo <- factor(genotab.m$haplo)
  genotab.m
}
```


```r
cnvfile1 <- "comb_FCs.CNVs.bed"
cnvs <- getGenoTab(cnvfile1)
cnvs$ypos <- as.numeric(cnvs$indv)
cnvs$indv <- factor(cnvs$indv,levels=tree_order)

cnvfile2 <- "comb_FCs.CNVs.coreGenome.bed"
cnvs_core <- getGenoTab(cnvfile2)
cnvs_core$ypos <- as.numeric(cnvs_core$indv)
cnvs_core$indv <- factor(cnvs_core$indv,levels=tree_order)
```


```r
cols=c("red","green")
names(cols)<-c("-","+")
colsc <- scale_color_manual(values = cols)

dkcols=c("dark red","dark green")
names(dkcols)<-c("-","+")
dkcolsc <- scale_color_manual(values = dkcols)
```



```r
#Clade1

#all vars
ggplot(cnvs,aes(xmin=ST,xmax=EN,ymin=indv,ymax=indv,colour=value)) + geom_rect(size=3) + 
  facet_grid(CHR ~ .) +colsc+ xlim(0,3.5e6) +ggtitle("CNVs, whole genome")
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-1.png) 

```r
ggplot(cnvs_core,aes(xmin=ST,xmax=EN,ymin=indv,ymax=indv,colour=value)) + geom_rect(size=3) + 
  facet_grid(CHR ~ .) +colsc+ xlim(0,3.5e6) +ggtitle("CNVs, core genome")
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-2.png) 

```r
ggplot(subset(cnvs,indv %in% c1),aes(xmin=ST,xmax=EN,ymin=indv,ymax=indv,colour=value)) + 
  geom_rect(size=2, colour="grey") + dkcolsc +
  geom_rect(data=subset(cnvs_core,indv %in% c1),size=3) + 
  facet_grid(CHR ~ .) +colsc+ xlim(0,3.5e6) +ggtitle("CNVs, clade1")
```

```
## Scale for 'colour' is already present. Adding another scale for 'colour', which will replace the existing scale.
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-3.png) 

```r
ggplot(subset(cnvs,indv %in% c2),aes(xmin=ST,xmax=EN,ymin=indv,ymax=indv,colour=value)) + 
  geom_rect(size=2, colour="grey") + dkcolsc +
  geom_rect(data=subset(cnvs_core,indv %in% c2),size=3) + 
  facet_grid(CHR ~ .) +colsc+ xlim(0,3.5e6) +ggtitle("CNVs, clade2")
```

```
## Scale for 'colour' is already present. Adding another scale for 'colour', which will replace the existing scale.
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-4.png) 

```r
ggplot(subset(cnvs,indv %in% c3),aes(xmin=ST,xmax=EN,ymin=indv,ymax=indv,colour=value)) + 
  geom_rect(size=2, alpha=0.5) + dkcolsc +
  geom_rect(data=subset(cnvs_core,indv %in% c3),size=3) + 
  facet_grid(CHR ~ .) +colsc+ xlim(0,3.5e6) +ggtitle("CNVs, clade3")
```

```
## Scale for 'colour' is already present. Adding another scale for 'colour', which will replace the existing scale.
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-5.png) 

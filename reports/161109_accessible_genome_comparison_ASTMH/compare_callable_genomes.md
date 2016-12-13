library(ggplot2)
library(knitr)
library(reshape2)
library(hexbin)

```r
opts_chunk$set(fig.width=12, fig.height=7, fig.path="figure/chrs-")
opts_chunk$set(dev=c('png','postscript'))
```


```r
calcols <- c("chr","st","en","NV","LDfails","LDpass","RSfails","RSpass","df","dn","Dpass","multipass")
callable1k_disco <- read.table("qDD2vr3D7.callable.blocks1k.tab.txt",stringsAsFactors = T,sep='\t',header=F,col.names = calcols)
callable1k_haplo <- read.table("fakeNGS_haplo_qDD2r3D7_rl100.callable.1k.txt",stringsAsFactors = T,sep='\t',header=F,col.names = calcols)

callable1k_disco <- callable1k_disco[,c(1,2,3,12)]
callable1k_haplo <- callable1k_haplo[,c(1,2,3,12)]
colnames(callable1k_disco)[[4]] <- "disco"
colnames(callable1k_haplo)[[4]] <- "haplo"
callable1k_haplo$haplo = callable1k_haplo$haplo=="True"
callable1k_disco$disco = callable1k_disco$disco=="True"

callable1k<-merge(callable1k_disco,callable1k_haplo)
# callable1k$accessible <- NULL
# callable1k[callable1k$disco=="True" & callable1k$haplo=="False","accessible"] <- "DISCOVAR only"
# callable1k[callable1k$disco=="False" & callable1k$haplo=="True","accessible"] <- "haplo only"
# callable1k[callable1k$disco=="False" & callable1k$haplo=="False","accessible"] <- "inacessible"

levels(callable1k$chr) <- c("MIT","01","02","03","04","05","06","07","08","09","10","11","12","13","14","API")
callable1k <- subset(callable1k,!chr %in% c("MIT","API"))
```




```r
calcols <- c("chr","st","en","NV","LDfails","LDpass","RSfails","RSpass","df","dn","Dpass","multipass")
callable5k_disco <- read.table("qDD2vr3D7.callable.blocks5k.tab.txt",stringsAsFactors = T,sep='\t',header=F,col.names = calcols)
callable5k_haplo <- read.table("fakeNGS_haplo_qDD2r3D7_rl100.callable.5k.txt",stringsAsFactors = T,sep='\t',header=F,col.names = calcols)

callable5k_disco <- callable5k_disco[,c(1,2,3,12)]
callable5k_haplo <- callable5k_haplo[,c(1,2,3,12)]

colnames(callable5k_disco)[[4]] <- "disco"
colnames(callable5k_haplo)[[4]] <- "haplo"
callable5k_haplo$haplo = callable5k_haplo$haplo=="True"
callable5k_disco$disco = callable5k_disco$disco=="True"

callable5k<-merge(callable5k_disco,callable5k_haplo)
# callable5k$accessible <- NULL
# callable5k[callable5k$disco=="True" & callable5k$haplo=="False","accessible"] <- "DISCOVAR only"
# callable5k[callable5k$disco=="False" & callable5k$haplo=="True","accessible"] <- "haplo only"
# callable5k[callable5k$disco=="False" & callable5k$haplo=="False","accessible"] <- "inacessible"

levels(callable5k$chr) <- c("MIT","01","02","03","04","05","06","07","08","09","10","11","12","13","14","API")
callable5k <- subset(callable5k,!chr %in% c("MIT","API"))
```



```r
telos <- readLines("List.subtelomeres.3D7.regions.txt")
telos <- t(as.data.frame(strsplit(telos,split = ':')))
telos <- as.data.frame(cbind(telos[,1],t(as.data.frame(strsplit(telos[,2],split = '-')))))
colnames(telos) <- c("chr","st","en")
telos$st <- as.numeric(as.character(telos$st))
telos$en <- as.numeric(as.character(telos$en))
levels(telos$chr) <- c("01","02","03","04","05","06","07","08","09","10","11","12","13","14")

telos <- merge(subset(telos,st < 1000),subset(telos,st > 1000),by="chr")
colnames(telos) <- c("chr","stL","enL","stR","enR")


#callable5k$chr <- as.character(callable5k$chr)
callable5k$telo=F
tmp <- merge(callable5k,telos,by="chr",all.x=T)

tmp[is.na(tmp$stL),c("stL","enL","stR","enR")] <- -1
tmp$telo[tmp$en < tmp$enL] = T # "L"
tmp$telo[tmp$st > tmp$stR] = T # "R"
callable5k <- tmp[,colnames(callable5k)]
rm(tmp)
```


```r
callable5k$accessible=F
callable5k$accessible[callable5k$disco & !callable5k$haplo] <- "DISCOVAR"
callable5k$accessible[callable5k$haplo & !callable5k$disco] <- "GATK"
callable5k$accessible[callable5k$haplo & callable5k$disco] <- NA

callable5k$accessible[!callable5k$disco & !callable5k$haplo] <- "INACCESSIBLE"

callable5k$accessible[!callable5k$disco & !callable5k$haplo & callable5k$telo] <- "INACCESSIBLE (TELOMERE)"

callable5k$accessible <- factor(callable5k$accessible,levels=c("DISCOVAR",
                                                               "GATK",
                                                               "INACCESSIBLE",
                                                               "INACCESSIBLE (TELOMERE)"
                                                               ))
```


```r
chrlen <- as.data.frame(aggregate(callable5k$en,list(callable5k$chr),FUN=max))
chrlen$st=1
colnames(chrlen) <- c("chr","st","en")
```



```r
blank_theme <-   theme(axis.line.y=element_blank(), axis.ticks.y=element_blank(),axis.text.y=element_blank(),
        panel.background=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="bottom")
ggplot(callable5k,aes(xmin=st,xmax=en,ymin=0,ymax=1)) + 
  ggtitle(paste("callable regions (5kb windows)")) +
  geom_rect(aes(fill=accessible)) + scale_fill_manual(values=c("blue","orange","grey","#888888"))+
  geom_rect(data=chrlen,fill=NA,colour="black") + 
  facet_grid(chr ~ .)   + theme(axis.text.y = element_blank(),axis.ticks.y=element_blank(),legend.position="bottom") +
  theme(panel.background = element_rect(fill='white'), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
```

![plot of chunk unnamed-chunk-6](figure/chrs-unnamed-chunk-6-1.png)

```r
ggplot(callable1k,aes(xmin=st,xmax=en,ymin=0,ymax=1)) + 
  ggtitle(paste("accessible regions (1kb windows)")) +
  geom_rect(aes(fill=accessible)) + scale_fill_manual(values=c("blue","orange","grey"))+
  geom_rect(data=chrlen, fill=NA,colour="black")+
  facet_grid(ch ~ .)  + blank_theme
```

```
## Error in layout_base(data, rows, drop = drop): At least one layer must contain all variables used for facetting
```

![plot of chunk unnamed-chunk-6](figure/chrs-unnamed-chunk-6-2.png)



```r
dust <- read.table("Pf3D7_v3.dust.bed.gz",col.names=c("chr","st","en"))
levels(dust$chr) <- c("MIT","01","02","03","04","05","06","07","08","09","10","11","12","13","14","API")
dust$len = dust$en-dust$st
callable5k$dustpc <- 0
#for (i in 1:20) {
for (i in 1:dim(callable5k)[1]) {
  mychr <- callable5k[i,1]
  myst <- callable5k[i,2]
  myen <- callable5k[i,3]
  dustset <- subset(dust,chr==mychr & en>=myst & st<myen)
  if (dim(dustset)[1]>0){
    dusted <- sum(dustset$len)
    olmin <- min(dustset$st)
    olmax <- max(dustset$en)
    if (olmin<st) {dusted<-dusted-(myst-olmin)}
    if (olmax>en) {dusted<-dusted-(olmax-myen)}
    dustpc = dusted/(en-st)
    callable5k[i,"dustpc"]<-dustpc  
    }
}
  
head(callable5k)
```

```
##   chr     st     en disco haplo telo              accessible    dustpc
## 1  01      1   5000 FALSE FALSE TRUE INACCESSIBLE (TELOMERE) 0.3674735
## 2  01 100001 105000  TRUE FALSE TRUE                DISCOVAR 0.5093019
## 3  01  10001  15000 FALSE FALSE TRUE INACCESSIBLE (TELOMERE) 0.2310462
## 4  01 105001 110000  TRUE  TRUE TRUE                    <NA> 0.4999000
## 5  01 110001 115000  TRUE  TRUE TRUE                    <NA> 0.2426485
## 6  01 115001 120000  TRUE  TRUE TRUE                    <NA> 0.5063013
```





```r
ggplot(callable5k,aes(xmin=st,xmax=en,ymin=0,ymax=1)) + 
  ggtitle(paste("callable regions (5kb windows)")) +
  geom_rect(aes(fill=accessible,colour=accessible)) + scale_fill_manual(values=c("blue","orange","grey","#888888"))+ scale_colour_manual(values=c("blue","orange","grey","#888888"))+
  geom_rect(data=chrlen,fill=NA,colour="black") + 
  facet_grid(chr ~ .)   + theme(axis.text.y = element_blank(),axis.ticks.y=element_blank(),legend.position="bottom") +
  theme(panel.background = element_rect(fill='white'), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
```

![plot of chunk unnamed-chunk-8](figure/chrs-unnamed-chunk-8-1.png)

```r
call5kM <- melt(callable5k,id.vars = c("chr","st","en","accessible","dustpc","telo"))

ggplot(subset(call5kM, value==T),aes(xmin=st,xmax=en,ymin=0,ymax=1)) + 
  geom_rect(aes(fill=variable)) + scale_fill_manual(values=c("blue","orange","grey","#888888"))+
  geom_rect(data=chrlen,fill=NA,colour="black") + 
  facet_grid(chr ~ variable)   + theme(axis.text.y = element_blank(),axis.ticks.y=element_blank(),legend.position="bottom") +
  theme(panel.background = element_rect(fill='white'), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
```

![plot of chunk unnamed-chunk-8](figure/chrs-unnamed-chunk-8-2.png)


```r
mychr="07"
callable5k$pos<-callable5k$st+(callable5k$en-callable5k$st)


accgr <- ggplot(subset(callable5k,chr==mychr),aes(xmin=st,xmax=en,ymin=0,ymax=1)) + 
  geom_rect(aes(fill=accessible,colour=accessible)) + scale_fill_manual(values=c("blue","orange","grey","#888888"))+ scale_colour_manual(values=c("blue","orange","grey","#888888"))+
  geom_rect(data=subset(chrlen,chr==mychr),fill=NA,colour="black") + 
  theme(axis.text.y = element_blank(),axis.ticks.y=element_blank(),legend.position="bottom") +
  theme(panel.background = element_rect(fill='white'), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

accgrsep <- ggplot(subset(call5kM, chr==mychr & value==T),aes(xmin=st,xmax=en,ymin=0,ymax=1)) + 
  geom_rect(aes(fill=variable)) + scale_fill_manual(values=c("blue","orange","grey","#888888"))+
  geom_rect(data=subset(chrlen,chr==mychr),fill=NA,colour="black") + 
  facet_grid(variable ~ .)   + theme(axis.text.y = element_blank(),axis.ticks.y=element_blank(),legend.position="none") +
  theme(panel.background = element_rect(fill='white'), panel.grid.major = element_blank(),strip.text.y = element_blank(), panel.grid.minor = element_blank())


dustgr <- ggplot(subset(callable5k,chr==mychr),aes(x=pos,y=dustpc)) + 
  geom_line()+xlab("")+ylab("repeat content")+xlim(0,1.5e6)+
  theme(axis.text.x=element_blank())

grid.arrange(dustgr, accgrsep,accgr, ncol=1,widths=c(1),heights=c(3,3,2))
```

![plot of chunk unnamed-chunk-9](figure/chrs-unnamed-chunk-9-1.png)
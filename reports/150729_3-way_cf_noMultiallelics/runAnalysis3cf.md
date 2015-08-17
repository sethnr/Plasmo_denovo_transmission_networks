library(ggplot2)
library(knitr)
library(reshape2)

```r
opts_chunk$set(fig.width=9, fig.height=6)
```


```r
depths <- read.table("cfNucmerDiscoHaplo_150729_leftAlignGATKSPLIT.DEPTH.txt",stringsAsFactors = F,sep='\t',header=T)
# colnames(values) <- c("chrom","start","end",
#                       "covered1","covered2","covered3",
#                       "d1",  "d2",  "d3",	
#                       "dcov1",  "dcov2",	"dcov3",
#                       "SNPS",	"Smatch0",	"Smatch1",	"Spriv",	
#                       "INDELS",	"Imatch",	"Imiss",	"Ipriv")
depths$block <- depths$end - ((depths$end-depths$start)/2)
depths$meandepth <- (depths$d2+depths$d3)/2
depths$block <- c(1:30) #nb, 10kb blocks
```


```r
nuccfs <- read.table("cfNucmerDiscoHaplo_150729_leftAlignGATKSPLIT.txt",sep="\t",stringsAsFactors = F)
colnames(nuccfs) <- c("file","chr","pos","type","concordance","matched","quality","length","alleles","complexity",
                      "STR","period","exponent","STRlength",
                      "STRcomplexity","Apc","Tpc","Cpc","Gpc" )
matchlevels <- c("",     "NUCMER","DISCO","HAPLO","NUCMER,DISCO","NUCMER,HAPLO","DISCO,HAPLO","NUCMER,DISCO,HAPLO")
colours <-     c("black","green", "red",  "blue", "yellow",      "cyan",        "magenta",     "white")
names(colours) <- matchlevels


nuccfs$matched <- factor(nuccfs$matched,levels = matchlevels)

nuccfs$polyA = (nuccfs$Apc==1 | nuccfs$Tpc==1)
nuccfs$TA = (nuccfs$Apc==0.5 & nuccfs$Tpc==0.5)
nuccfs$STRtype=""
nuccfs$STRtype <- factor(nuccfs$STRtype,levels=c("","STR","TA","polyA"))
nuccfs[nuccfs$STR != "NULL",]$STRtype = "STR"
nuccfs[nuccfs$TA,]$STRtype = "TA"
nuccfs[nuccfs$polyA,]$STRtype = "polyA"

nuccfs$quality <- as.numeric(nuccfs$quality)
nuccfs <- nuccfs[which(nuccfs$concordance!="MISMATCH"),]
```

#method concordance
##left align

```r
ctable <- as.data.frame(table(subset(nuccfs,concordance != "MATCH0",c(matched,type,STRtype))))
ctable$x <- c(0,1,2,3,1,1,2,4)
ctable$y <- c(0,1,2,3,2,3,3,4)
ctable$x <- factor(ctable$x,labels=c("none","NUCMER","DISCO","HAPLO","ALL"))
ctable$y <- factor(ctable$y,labels=c("none","NUCMER","DISCO","HAPLO","ALL"))
ggplot(subset(ctable,x!="none"),aes(x,y,fill=Freq)) + 
  geom_tile() + scale_fill_gradient(low="white",high="red") +
  stat_bin2d(geom="text", aes(label=Freq)) +
  facet_grid(type ~ STRtype) + ggtitle("concordance, aligned indels")
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-1.png) 

#call concordance vs posn, SNPs

```r
posConcL <- ggplot(subset(nuccfs,type=="S"),aes(x=pos, group=matched, fill=matched)) + geom_histogram(binwidth=10000)
posConcL + ggtitle("concordance x pos, aligned SNPs") + scale_fill_manual(values = colours)
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-1.png) 
#call concordance vs posn, INDELS

```r
posConcL <- ggplot(subset(nuccfs,type=="I"),aes(x=pos, group=matched, fill=matched)) + geom_histogram(binwidth=10000)
posConcL + ggtitle("concordance x pos, aligned INDELS") + scale_fill_manual(values = colours)
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-1.png) 
##STR type (A/TA/complex) vs concordance, INDELS

```r
posConcSTRL <- ggplot(subset(nuccfs,type=="I"),aes(x=pos, group=matched, fill=matched)) + geom_histogram(binwidth=10000)  + facet_grid(STRtype ~ .,scales="free_y") + ggtitle("STR type x pos, aligned INDELS") + scale_fill_manual(values = colours)
posConcSTRL
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6-1.png) 
library(ggplot2)
library(knitr)

```r
#opts_chunk$set(fig.width=9, fig.height=6)
opts_chunk$set(fig.width=9, fig.height=4.5)
```

#depth vs snp rates vs call rates

```r
values <- read.table("depth_v_het_10kb_150629.txt",stringsAsFactors = F)
colnames(values) <- c("chrom","start","end",
                      "covered1","covered2","covered3",
                      "d1",  "d2",	"d3",	
                      "dcov1",  "dcov2",	"dcov3",
                      "SNPS",	"Smatch",	"Smiss",	"Spriv",	"INDELS",	"Imatch",	"Imiss",	"Ipriv")
values$block <- values$end - ((values$end-values$start)/2)
values$meandepth <- (values$d2+values$d3)/2
```

##alignment depth v 3D7

```r
depths <- values[,c("chrom","block","d1","d2","d3")]
colnames(depths) <- c("chrom","block","3D7","DD2","IT")
depths <- melt(depths,id.vars =c("chrom","block"))

depthplot <- ggplot(depths,aes(x=block,y=value,colour=variable,group=variable)) + geom_line()+
  ylab("depth")+xlab("pos")
depthplot
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2-1.png) 

##INDEL counts & concordance w. nucmer

```r
indelcounts <- values[,c("chrom","block","INDELS",  "Imatch",	"Imiss",	"Ipriv")]
indelcounts <- melt(indelcounts,id.vars =c("chrom","block","INDELS"))
indelcounts <- droplevels(indelcounts)
ggplot(indelcounts,aes(x=block,y=value*INDELS,group=block,fill=variable)) + geom_bar(stat="identity")
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-1.png) 

```r
indelplot <- ggplot(indelcounts,aes(x=block,y=value,group=block,fill=variable)) + geom_bar(stat="identity") + ylab("indel concordance")
indelplot
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-2.png) 

##SNP counts & concordance w. nucmer

```r
snpcounts <- values[,c("chrom","block","SNPS",  "Smatch",  "Smiss",	"Spriv")]
snpcounts <- melt(snpcounts,id.vars =c("chrom","block","SNPS"))
snpcounts <- droplevels(snpcounts)
ggplot(snpcounts,aes(x=block,y=value*SNPS,group=block,fill=variable)) + geom_bar(stat="identity")
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-1.png) 

```r
snpplot <- ggplot(snpcounts,aes(x=block,y=value,group=block,fill=variable)) + geom_bar(stat="identity") + ylab("snp concordance")
snpplot
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-2.png) 



```r
#tracks(depthplot,indelplot,snpplot)
```

##matched indel proportions, SNPs v INDELs

```r
varcf <- values[,c("chrom","block","meandepth","SNPS",  "INDELS","Smatch",  "Imatch")]
ggplot(varcf,aes(x=Smatch,y=Imatch,colour=meandepth)) + geom_point()
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6-1.png) 

```r
#ggplot(varcf,aes(x=Smatch*SNPS,y=Imatch*INDELS,colour=meandepth)) + geom_point()
```

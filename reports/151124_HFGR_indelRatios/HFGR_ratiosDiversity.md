library(ggplot2)
library(reshape2)
library(knitr)
library(zoo)



```r
opts_chunk$set(fig.width=14, fig.height=8)
```



```r
window=10
modeI <- function(x) {as.numeric(names(sort(-table(as.numeric(x))))[1])}
mode <- function(x) {names(sort(-table(as.numeric(x))))[1]}


HF2alleles <- read.table("HFGR_all.incDD2.PASS.HFGR.LMRG.HAP.DISCORD_HF2.alleles.tab.txt",sep="\t",header=T)
HF2alleles$cumpos <- HF2alleles$pos + chrLens[as.character(HF2alleles$chr),"increment"]

HF2ratios <- data.frame("pos"=rollapply(HF2alleles$cumpos,FUN=mean,width=window),
  "ratio"=rollapply(HF2alleles$type,FUN=function(x) {sum(x=="SNP")/sum(x=="INDEL")},width=window),
  "chrNo"=rollapply(HF2alleles$chr,FUN=function(x) {length(unique(x))},width=window),
  "chr"=rollapply(HF2alleles$chr,FUN=function(x) {paste(as.character(unique(x)),sep=",",collapse=",")},width=window),
  "chrCol"=rollapply(HF2alleles$chr,FUN=function(x) {modeI(x) %% 2},width=window),
  "dataset"="HF2"
)

HF3alleles <- read.table("HFGR_all.incDD2.PASS.HFGR.LMRG.HAP.DISCORD_HF3.alleles.tab.txt",sep="\t",header=T)
HF3alleles$cumpos <- HF3alleles$pos + chrLens[as.character(HF3alleles$chr),"increment"]

HF3ratios <- data.frame("pos"=rollapply(HF3alleles$cumpos,FUN=mean,width=window),
  "ratio"=rollapply(HF3alleles$type,FUN=function(x) {sum(x=="SNP")/sum(x=="INDEL")},width=window),
  "chrNo"=rollapply(HF3alleles$chr,FUN=function(x) {length(unique(x))},width=window),
  "chr"=rollapply(HF3alleles$chr,FUN=function(x) {paste(as.character(unique(x)),sep=",",collapse=",")},width=window),
  "chrCol"=rollapply(HF3alleles$chr,FUN=function(x) {modeI(x) %% 2},width=window),
  "dataset"="HF3"
)
```


#HF2 and HF3 clades

```r
ratios <-rbind(allratios,HF2ratios,HF3ratios)
#ratios <- subset(ratios,chrNo==1)
ggplot(ratios,aes(x=pos,y=ratio,colour=chrCol)) + geom_point() + 
  geom_hline(aes(yintercept=2)) + geom_vline(aes(xintercept=5e5+chrLens[as.character("Pf3D7_12_v3"),"increment"])) + 
  scale_y_log10() + facet_grid(dataset ~ .) + xlim(0,24e6)
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2-1.png) 

#HF3 clades only

```r
ratios <-rbind(clade1ratios,clade2ratios,clade3ratios)
ratios <- subset(ratios,chrNo==1)
ggplot(ratios,aes(x=pos,y=ratio,colour=chrCol)) + geom_point() + geom_hline(aes(yintercept=2)) + scale_y_log10() + facet_grid(dataset ~ .) + xlim(0,24e6)
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-1.png) 

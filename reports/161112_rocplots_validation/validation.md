library(ggplot2)
library(hexbin)
library(reshape2)
library(knitr)

library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)

```r
opts_chunk$set(fig.width=12, fig.height=8,dev=c('png','postscript'),warning=F)
```


```r
discoFilt <- read.table("3D7DD2.REFCALL.PASS.DD2CONC.PFX.LMRG.HAP.Dd2-2D4_FILT_EVAL/weighted_roc.tsv.gz",header=F,skip = 3,stringsAsFactors =F,sep="\t",col.names=c("score","true_positives","false_positives","false_negatives","precision","sensitivity","f_measure"))
discoFilt$caller="disco"
discoFilt$region="GENOME"

discoFiltSNPs <- read.table("3D7DD2.REFCALL.PASS.DD2CONC.PFX.LMRG.HAP.Dd2-2D4_FILT_EVAL/snp_roc.tsv.gz",header=F,skip = 3,stringsAsFactors =F,sep="\t",col.names = c("score","true_positives","false_positives"))
discoFiltINDELs <- read.table("3D7DD2.REFCALL.PASS.DD2CONC.PFX.LMRG.HAP.Dd2-2D4_FILT_EVAL/non_snp_roc.tsv.gz",header=F,skip = 3,stringsAsFactors =F,sep="\t",col.names = c("score","true_positives","false_positives"))
discoFiltSNPs$type="SNP"
discoFiltINDELs$type="INDEL"
discoFiltType <- rbind(discoFiltSNPs,discoFiltINDELs)
discoFiltType$caller="disco"
discoFiltType$region="GENOME"
posns <- intersect(discoFiltSNPs$score,discoFiltINDELs$score)
discoFiltType <- subset(discoFiltType,score %in% posns)


discoCALLFilt <- read.table("3D7DD2.REFCALL.PASS.DD2CONC.PFX.LMRG.HAP.CALLBOTH.Dd2-2D4_FILT_EVAL/weighted_roc.tsv.gz",header=F,skip = 3,stringsAsFactors =F,sep="\t",col.names=c("score","true_positives","false_positives","false_negatives","precision","sensitivity","f_measure"))
discoCALLFilt$caller="disco"
discoCALLFilt$region="CALLABLE"

discoCALLFiltSNPs <- read.table("3D7DD2.REFCALL.PASS.DD2CONC.PFX.LMRG.HAP.CALLBOTH.Dd2-2D4_FILT_EVAL/snp_roc.tsv.gz",header=F,skip = 3,stringsAsFactors =F,sep="\t",col.names = c("score","true_positives","false_positives"))
discoCALLFiltINDELs <- read.table("3D7DD2.REFCALL.PASS.DD2CONC.PFX.LMRG.HAP.CALLBOTH.Dd2-2D4_FILT_EVAL/non_snp_roc.tsv.gz",header=F,skip = 3,stringsAsFactors =F,sep="\t",col.names = c("score","true_positives","false_positives"))
discoCALLFiltSNPs$type="SNP"
discoCALLFiltINDELs$type="INDEL"
discoCALLFiltType <- rbind(discoCALLFiltSNPs,discoCALLFiltINDELs)
discoCALLFiltType$caller="disco"
discoCALLFiltType$region="CALLABLE"
posns <- intersect(discoCALLFiltSNPs$score,discoCALLFiltINDELs$score)
discoCALLFiltType <- subset(discoCALLFiltType,score %in% posns)


head(discoFilt); head(discoCALLFilt)
```

```
##     score true_positives false_positives false_negatives precision
## 1 5539.62              1               0           19976         1
## 2 4964.62              2               0           19975         1
## 3 4904.62              3               0           19974         1
## 4 4788.50              4               0           19973         1
## 5 4701.25              5               0           19972         1
## 6 4553.38              6               0           19971         1
##   sensitivity f_measure caller region
## 1       1e-04     1e-04  disco GENOME
## 2       1e-04     2e-04  disco GENOME
## 3       2e-04     3e-04  disco GENOME
## 4       2e-04     4e-04  disco GENOME
## 5       3e-04     5e-04  disco GENOME
## 6       3e-04     6e-04  disco GENOME
```

```
##     score true_positives false_positives false_negatives precision
## 1 4964.62              1               0           11864         1
## 2 4553.38              2               0           11863         1
## 3 3912.62              3               0           11862         1
## 4 3860.75              4               0           11861         1
## 5 3721.12              5               0           11860         1
## 6 3685.62              6               0           11859         1
##   sensitivity f_measure caller   region
## 1       1e-04     2e-04  disco CALLABLE
## 2       2e-04     3e-04  disco CALLABLE
## 3       3e-04     5e-04  disco CALLABLE
## 4       3e-04     7e-04  disco CALLABLE
## 5       4e-04     8e-04  disco CALLABLE
## 6       5e-04     1e-03  disco CALLABLE
```

```r
head(discoFiltType); head(discoCALLFiltType)
```

```
##       score true_positives false_positives type caller region
## 543 2858.50            568               5  SNP  disco GENOME
## 602 2810.88            630               6  SNP  disco GENOME
## 655 2767.62            687               7  SNP  disco GENOME
## 685 2743.00            722               9  SNP  disco GENOME
## 761 2694.00            809              12  SNP  disco GENOME
## 796 2671.25            850              12  SNP  disco GENOME
```

```
##       score true_positives false_positives type caller   region
## 489 2810.88            509               4  SNP  disco CALLABLE
## 532 2767.62            554               5  SNP  disco CALLABLE
## 619 2694.00            649              10  SNP  disco CALLABLE
## 681 2627.25            726              12  SNP  disco CALLABLE
## 690 2617.00            735              12  SNP  disco CALLABLE
## 703 2608.12            749              12  SNP  disco CALLABLE
```

```r
# discoFilt <- read.table("disco_3D7_0901_weighted_roc.tsv.gz",header=F,skip = 3,stringsAsFactors =F,sep="\t")
# colnames(discoFilt) <-c("score","true_positives","false_positives","false_negatives","precision","sensitivity","f_measure")
# discoFilt$sample="0901"
# discoFilt$caller="disco"
# discoFilt$filter=T
```


```r
gatkFilt <- read.table("3D7DD2_300100.VQSR.Pass90.RENAME.Dd2-2D4_FILT_EVAL/weighted_roc.tsv.gz",header=F,skip = 3,stringsAsFactors =F,sep="\t",col.names=c("score","true_positives","false_positives","false_negatives","precision","sensitivity","f_measure"))
gatkFilt$caller="gatk"
gatkFilt$region="GENOME"

gatkFiltSNPs <- read.table("3D7DD2_300100.VQSR.Pass90.RENAME.Dd2-2D4_FILT_EVAL/snp_roc.tsv.gz",header=F,skip = 3,stringsAsFactors =F,sep="\t",col.names = c("score","true_positives","false_positives"))
gatkFiltINDELs <- read.table("3D7DD2_300100.VQSR.Pass90.RENAME.Dd2-2D4_FILT_EVAL/non_snp_roc.tsv.gz",header=F,skip = 3,stringsAsFactors =F,sep="\t",col.names = c("score","true_positives","false_positives"))
gatkFiltSNPs$type="SNP"
gatkFiltINDELs$type="INDEL"
gatkFiltType <- rbind(gatkFiltSNPs,gatkFiltINDELs)
gatkFiltType$caller="gatk"
gatkFiltType$region="GENOME"
posns <- intersect(gatkFiltSNPs$score,gatkFiltINDELs$score)
gatkFiltType <- subset(gatkFiltType,score %in% posns)


gatkCALLFilt <- read.table("3D7DD2_300100.VQSR.Pass90.RENAME.CALLBOTH.Dd2-2D4_FILT_EVAL/weighted_roc.tsv.gz",header=F,skip = 3,stringsAsFactors =F,sep="\t",col.names=c("score","true_positives","false_positives","false_negatives","precision","sensitivity","f_measure"))
gatkCALLFilt$caller="gatk"
gatkCALLFilt$region="CALLABLE"

gatkCALLFiltSNPs <- read.table("3D7DD2_300100.VQSR.Pass90.RENAME.CALLBOTH.Dd2-2D4_FILT_EVAL/snp_roc.tsv.gz",header=F,skip = 3,stringsAsFactors =F,sep="\t",col.names = c("score","true_positives","false_positives"))
gatkCALLFiltINDELs <- read.table("3D7DD2_300100.VQSR.Pass90.RENAME.CALLBOTH.Dd2-2D4_FILT_EVAL/non_snp_roc.tsv.gz",header=F,skip = 3,stringsAsFactors =F,sep="\t",col.names = c("score","true_positives","false_positives"))
gatkCALLFiltSNPs$type="SNP"
gatkCALLFiltINDELs$type="INDEL"
gatkCALLFiltType <- rbind(gatkCALLFiltSNPs,gatkCALLFiltINDELs)
gatkCALLFiltType$caller="gatk"
gatkCALLFiltType$region="CALLABLE"
posns <- intersect(gatkCALLFiltSNPs$score,gatkCALLFiltINDELs$score)
gatkCALLFiltType <- subset(gatkCALLFiltType,score %in% posns)
```


```r
rocCf <- rbind(discoCALLFiltType,discoFiltType,gatkCALLFiltType,gatkFiltType)
```


```r
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
```



```r
#ROC PLOT
ggplot(discoCALLFiltType,aes(x=false_positives,y=true_positives,colour=type)) + 
  geom_line() + facet_grid(. ~ caller) + theme(legend.position="bottom")
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-1.png)

```r
ggplot(discoFiltType,aes(x=false_positives,y=true_positives,colour=type)) + 
  geom_point() + facet_grid(. ~ caller) + theme(legend.position="bottom")
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-2.png)

```r
ggplot(discoFiltType,aes(x=false_positives,y=true_positives,colour=type)) + 
  geom_line() + facet_grid(. ~ caller) + theme(legend.position="bottom")
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-3.png)

```r
ggplot(rocCf,aes(x=false_positives,y=true_positives,colour=type,linetype=region)) + 
  geom_line() + facet_grid(. ~ caller) + theme(legend.position="bottom")
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-4.png)

```r
ggplot(rocCf,aes(x=false_positives,y=true_positives,colour=type,linetype=region)) + 
  geom_line() + scale_x_log10() + facet_grid(. ~ caller) + theme(legend.position="bottom")
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-5.png)



```r
callboth <- read.table("Thies_all_manual.PASS.Cls.miss0.5.LMRG.HAP.CALLBOTH.vcf.dist.tab.txt",sep="\t")
callboth$from=rownames(callboth)
callboth <- melt(callboth,variable.name = "to",value.name="both")
```

```
## Using from as id variables
```

```r
callDisco <- read.table("Thies_all_manual.PASS.Cls.miss0.5.LMRG.HAP.CALLDISCOONLY.vcf.dist.tab.txt",sep="\t")
callDisco$from=rownames(callDisco)
callDisco <- melt(callDisco,variable.name = "to",value.name="disco")
```

```
## Using from as id variables
```

```r
callsnp <- read.table("Thies_all_manual.PASS.Cls.miss0.5.LMRG.HAP.SNP.recode.vcf.dist.tab.txt",sep="\t")
callsnp $from=rownames(callsnp )
callsnp  <- melt(callsnp ,variable.name = "to",value.name="snp")
```

```
## Using from as id variables
```

```r
callindel <- read.table("Thies_all_manual.PASS.Cls.miss0.5.LMRG.HAP.INDEL.recode.vcf.dist.tab.txt",sep="\t")
callindel $from=rownames(callindel )
callindel  <- melt(callindel ,variable.name = "to",value.name="indel")
```

```
## Using from as id variables
```

```r
callCf <- merge(merge(callboth,callDisco),merge(callsnp,callindel))
callCf <- callCf[!is.na(callCf$snp),]
callCf <- callCf[callCf$from != callCf$to,]

callCf <- callCf[callCf$snp <= 1000,]

callbothGATK <- read.table("thies_300100_haplo.CALLBOTH.RENAME.dist.tab.txt",sep="\t")
callbothGATK$from=rownames(callbothGATK)
callbothGATK <- melt(callbothGATK,variable.name = "to",value.name="both")
```

```
## Using from as id variables
```

```r
callhaploGATK <- read.table("thies_300100_haplo.CALLHAPLO.RENAME.dist.tab.txt",sep="\t")
callhaploGATK$from=rownames(callhaploGATK)
callhaploGATK <- melt(callhaploGATK,variable.name = "to",value.name="haplo")
```

```
## Using from as id variables
```

```r
calldiscoGATK <- read.table("thies_300100_haplo.CALLDISCO.RENAME.dist.tab.txt",sep="\t")
calldiscoGATK$from=rownames(calldiscoGATK)
calldiscoGATK <- melt(calldiscoGATK,variable.name = "to",value.name="disco")
```

```
## Using from as id variables
```

```r
callsnpGATK <- read.table("thies_300100_haplo.CALLHAPLO.RENAME.SNP.dist.tab.txt",sep="\t")
callsnpGATK$from=rownames(callsnpGATK)
callsnpGATK <- melt(callsnpGATK,variable.name = "to",value.name="snp")
```

```
## Using from as id variables
```

```r
callindelGATK <- read.table("thies_300100_haplo.CALLHAPLO.RENAME.INDEL.dist.tab.txt",sep="\t")
callindelGATK$from=rownames(callindelGATK)
callindelGATK <- melt(callindelGATK,variable.name = "to",value.name="indel")
```

```
## Using from as id variables
```

```r
callCfGATK <- merge(merge(callsnpGATK,callindelGATK),merge(callbothGATK,callhaploGATK))
callCfGATK <- merge(callCfGATK,calldiscoGATK)
callCfGATK <- callCfGATK[!is.na(callCfGATK$both),]
callCfGATK <- callCfGATK[callCfGATK$from != callCfGATK$to,]
head(callCfGATK)
```

```
##        from       to  snp indel  both haplo disco
## 6  Th061.13 Th095.13  435   740  1048  1175  1425
## 19 Th068.12 Th061.13  384   739  1012  1123  1356
## 24 Th068.12 Th095.13  404   782  1054  1186  1454
## 37 Th074.13 Th061.13 9698 12608 20163 22306 26425
## 38 Th074.13 Th068.12 9644 12591 20069 22235 26353
## 41 Th074.13 Th092.13 9790 12668 20335 22458 26871
```

```r
callCfGATK <- callCfGATK[callCfGATK$both > 0,]

callCfGATK <- callCfGATK[callCfGATK$both <= 5000,]

#remove core calls from disco/haplo regions
# discovar already called on disco only
#callCf$disco <- callCf$disco-callCf$both

callCfGATK$disco <- callCfGATK$disco-callCfGATK$both
callCfGATK$haplo <- callCfGATK$haplo-callCfGATK$both
```


```r
cor(callCf$snp,callCf$indel)
```

```
## [1] 0.9519626
```

```r
ggplot(callCf,aes(y=snp,x=indel)) + geom_point() +  
  geom_label(y=100,x=150,label=paste("r = ",round(cor(callCf$snp,callCf$indel),2)),size=5) +
  ggtitle("Discovar calls - SNP:INDEL correlation")
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-1.png)

```r
cor(callCf$both,callCf$disco)
```

```
## [1] 0.9062925
```

```r
ggplot(callCf,aes(y=both,x=disco)) + geom_point() +  
  geom_label(y=100,x=150,label=paste("r = ",round(cor(callCf$both,callCf$disco),2)),size=5) +
  ggtitle("Discovar calls - accessible:inaccessible correlation")
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-2.png)

```r
#ggplot(callCf,aes(x=both,y=haplo)) + geom_point(aes(y=disco),colour="red")

cor(callCfGATK$both,callCfGATK$disco)
```

```
## [1] 0.9854283
```

```r
cor(callCfGATK$both,callCfGATK$haplo)
```

```
## [1] 0.9700306
```

```r
ggplot(callCfGATK,aes(y=both,x=haplo)) + 
  geom_point(aes(x=haplo),colour="blue") + geom_point(aes(x=disco),colour="red") +
  geom_label(y=1500,x=375,label=paste("r = ",round(cor(callCfGATK$both,callCfGATK$haplo),2)),size=5,colour="blue") + geom_label(y=1000,x=750,label=paste("r = ",round(cor(callCfGATK$both,callCfGATK$disco),2)),size=5,colour="red") +
  ggtitle("GATK calls - accessible:inaccessible correlation")
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-3.png)

```r
cor(callCfGATK$snp,callCfGATK$indel)
```

```
## [1] 0.8925991
```

```r
ggplot(callCfGATK,aes(y=snp,x=indel)) + geom_point() + 
  geom_label(x=1500,y=300,label=paste("r = ",round(cor(callCfGATK$snp,callCfGATK$indel),2)),size=5) +
  ggtitle("GATK calls - SNP:INDEL correlation")
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-4.png)

```r
callCfGATK$caller<-"GATK"
callCf$caller<-"DISCOVAR"
callCfAll <- rbind(callCfGATK[,c("from","to","caller","snp","indel")],callCf[,c("from","to","caller","snp","indel")])
ggplot(callCfAll,aes(y=snp,x=indel)) + geom_point() + geom_smooth(method = lm,se = F,linetype=2) +
  facet_grid(. ~ caller,scale="free_x")+
  ggtitle("SNP:INDEL correlation")+
  theme(legend.position="bottom")
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-5.png)


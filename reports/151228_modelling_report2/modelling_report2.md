library(ggplot2)
library(reshape2)
library(knitr)
library(splines)

```r
opts_chunk$set(fig.width=12, fig.height=6)
```


```r
onemth <- read.table("all_AFs_28days.txt",sep="\t",row.names=1)
sixmth <- read.table("all_AFs_180days.txt",sep="\t",row.names=1)

afcfs <- cbind(t(sixmth),t(onemth))

afcfs_m <- melt(afcfs)

keys <- t(as.data.frame(strsplit(as.character(afcfs_m$Var2),"_")))
afcfs_m <- cbind(keys, afcfs_m[,2:3])
colnames(afcfs_m) <- c("run","days","id","MAF")

cols<-c("90"="blue","14"="red")

afcfs_m <- afcfs_m[afcfs_m$MAF!=0,]

key <- unique(afcfs_m[,3:1])
key$days <- as.character(key$days)
key$colour=cols[key$days]

idcols <- key$colour
names(idcols) <- key$id
daycols <- key$colour
names(daycols) <- key$day
```


```r
ggplot(afcfs_m,aes(x=MAF,group=id,colour=days)) + geom_density(adjust=0.5) + xlim(0,0.1) +  ggtitle("sfs short/ long term infections") + scale_colour_manual(values=daycols)
```

```
## Warning in loop_apply(n, do.ply): Removed 8 rows containing non-finite
## values (stat_density).
```

```
## Warning in loop_apply(n, do.ply): Removed 18 rows containing non-finite
## values (stat_density).
```

```
## Warning in loop_apply(n, do.ply): Removed 3 rows containing non-finite
## values (stat_density).
```

```
## Warning in loop_apply(n, do.ply): Removed 7 rows containing non-finite
## values (stat_density).
```

```
## Warning in loop_apply(n, do.ply): Removed 22 rows containing non-finite
## values (stat_density).
```

```
## Warning in loop_apply(n, do.ply): Removed 8 rows containing non-finite
## values (stat_density).
```

```
## Warning in loop_apply(n, do.ply): Removed 17 rows containing non-finite
## values (stat_density).
```

```
## Warning in loop_apply(n, do.ply): Removed 9 rows containing non-finite
## values (stat_density).
```

```
## Warning in loop_apply(n, do.ply): Removed 11 rows containing non-finite
## values (stat_density).
```

```
## Warning in loop_apply(n, do.ply): Removed 1 rows containing non-finite
## values (stat_density).
```

```
## Warning in loop_apply(n, do.ply): Removed 9 rows containing non-finite
## values (stat_density).
```

```
## Warning in loop_apply(n, do.ply): Removed 6 rows containing non-finite
## values (stat_density).
```

```
## Warning in loop_apply(n, do.ply): Removed 12 rows containing non-finite
## values (stat_density).
```

```
## Warning in loop_apply(n, do.ply): Removed 1 rows containing non-finite
## values (stat_density).
```

```
## Warning in loop_apply(n, do.ply): Removed 1 rows containing non-finite
## values (stat_density).
```

```
## Warning in loop_apply(n, do.ply): Removed 2 rows containing non-finite
## values (stat_density).
```

```
## Warning in loop_apply(n, do.ply): Removed 1 rows containing non-finite
## values (stat_density).
```

```
## Warning in loop_apply(n, do.ply): Removed 1 rows containing non-finite
## values (stat_density).
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2-1.png) 

```r
ggplot(subset(afcfs_m,MAF>0.01),aes(x=MAF,group=id,colour=days)) + geom_density(adjust=0.8) + xlim(0,0.1)+scale_colour_manual(values=daycols) + ggtitle("sfs  short/ long term infections (af > 1pc)")
```

```
## Warning in loop_apply(n, do.ply): Removed 8 rows containing non-finite
## values (stat_density).
```

```
## Warning in loop_apply(n, do.ply): Removed 18 rows containing non-finite
## values (stat_density).
```

```
## Warning in loop_apply(n, do.ply): Removed 3 rows containing non-finite
## values (stat_density).
```

```
## Warning in loop_apply(n, do.ply): Removed 7 rows containing non-finite
## values (stat_density).
```

```
## Warning in loop_apply(n, do.ply): Removed 22 rows containing non-finite
## values (stat_density).
```

```
## Warning in loop_apply(n, do.ply): Removed 8 rows containing non-finite
## values (stat_density).
```

```
## Warning in loop_apply(n, do.ply): Removed 17 rows containing non-finite
## values (stat_density).
```

```
## Warning in loop_apply(n, do.ply): Removed 9 rows containing non-finite
## values (stat_density).
```

```
## Warning in loop_apply(n, do.ply): Removed 11 rows containing non-finite
## values (stat_density).
```

```
## Warning in loop_apply(n, do.ply): Removed 1 rows containing non-finite
## values (stat_density).
```

```
## Warning in loop_apply(n, do.ply): Removed 9 rows containing non-finite
## values (stat_density).
```

```
## Warning in loop_apply(n, do.ply): Removed 6 rows containing non-finite
## values (stat_density).
```

```
## Warning in loop_apply(n, do.ply): Removed 12 rows containing non-finite
## values (stat_density).
```

```
## Warning in loop_apply(n, do.ply): Removed 1 rows containing non-finite
## values (stat_density).
```

```
## Warning in loop_apply(n, do.ply): Removed 1 rows containing non-finite
## values (stat_density).
```

```
## Warning in loop_apply(n, do.ply): Removed 2 rows containing non-finite
## values (stat_density).
```

```
## Warning in loop_apply(n, do.ply): Removed 1 rows containing non-finite
## values (stat_density).
```

```
## Warning in loop_apply(n, do.ply): Removed 1 rows containing non-finite
## values (stat_density).
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2-2.png) 

```r
ggplot(subset(afcfs_m,MAF>0.01),aes(x=MAF,group=days,colour=days)) + geom_density(adjust=1.5) + xlim(0,0.5)+scale_colour_manual(values=daycols) + ggtitle("sfs  short/ long term infections (af > 1pc)")
```

```
## Warning in loop_apply(n, do.ply): Removed 2 rows containing non-finite
## values (stat_density).
```

```
## Warning in loop_apply(n, do.ply): Removed 22 rows containing non-finite
## values (stat_density).
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2-3.png) 



```r
histo <- read.table("all_histograms.txt",header=F,sep="\t")
edges <- seq(0,0.25,length.out=20)
colnames(histo) <- c("run","gen","N","maxN",edges,1)

histo$pc2 <- apply(histo[,7:25],1,FUN=sum)
histo$pc5 <- apply(histo[,9:25],1,FUN=sum)
histo$pc10 <- apply(histo[,13:25],1,FUN=sum)

ggplot(histo,aes(x=gen,y=pc2,colour=run)) + geom_point()
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-1.png) 

```r
detects <- aggregate(histo[,c("pc2","pc5","pc10")],list("gen"=histo$gen),FUN=mean)
detects <- melt(detects,id.vars=c("gen"))
colnames(detects) <- c("gen","MAF","mean")
levels(detects$MAF) <- c(">0.02",">0.05",">0.1")
detects$days <- detects$gen*2

#LOESS
ggplot(detects,aes(x=days,y=mean,colour=MAF))+ geom_point() +
  geom_smooth(method = "loess" ) + 
   geom_vline(xintercept=28,linetype=2) + geom_vline(xintercept=180,linetype=2)+
  ggtitle("LOESS")
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-2.png) 

```r
#GLM
ggplot(detects,aes(x=days,y=mean,colour=MAF))+ geom_point() +
#  geom_smooth(method = "lm", formula = y ~ ns(x,3)) + 
  geom_smooth(method = "glm", formula = y ~ ns(x,4) ) + 
   geom_vline(xintercept=28,linetype=2) + geom_vline(xintercept=180,linetype=2)+
  ggtitle("GLM")
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-3.png) 

```r
#LOESS, lim 120 days
ggplot(detects,aes(x=days,y=mean,colour=MAF))+ geom_point() +
  geom_smooth(method = "loess") + xlim(0,130) +
   geom_vline(xintercept=28,linetype=2) + geom_vline(xintercept=120,linetype=2)+
  ggtitle("LOESS, lim 130, samp 120")
```

```
## Warning in loop_apply(n, do.ply): Removed 26 rows containing missing
## values (stat_smooth).
```

```
## Warning in loop_apply(n, do.ply): Removed 26 rows containing missing
## values (stat_smooth).
```

```
## Warning in loop_apply(n, do.ply): Removed 26 rows containing missing
## values (stat_smooth).
```

```
## Warning in loop_apply(n, do.ply): Removed 78 rows containing missing
## values (geom_point).
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-4.png) 
library(ggplot2)
library(reshape2)
library(knitr)

```r
opts_chunk$set(fig.width=12, fig.height=6)
```


```r
afs_small <- data.frame(t(read.table("qstages.C1.5m.L0.01mb.G36.mtp.AFs.txt",sep="\t",row.names=1)))
colnames(afs_small)<-c("mth","halfyr","yr")
afs_large <- data.frame(t(read.table("qstages.C100.0m.L0.01mb.G36.mtp.AFs.txt",sep="\t",row.names=1)))
colnames(afs_large)<-c("mth","halfyr","yr")
# plot(density(afs_large$mth),col="red") + lines(density(afs_small$yr),col="blue")
# 
# 
# plot(density(subset(afs_large,mth>0.01)$mth),col="red",xlim=c(0,0.1)) + 
#   lines(density(subset(afs_small,halfyr>0.01)$halfyr),col="blue") +
#   lines(density(subset(afs_small,yr>0.01)$yr),col="dark blue")
# 

afcfs <- cbind(afs_large$mth,afs_small$halfyr)
colnames(afcfs) <- c("short_term","long_term")
afcfs_m <- melt(afcfs)

cols<-c("long_term"="blue","short_term"="red")

ggplot(afcfs_m,aes(x=value,group=Var2,colour=Var2)) + geom_density(adjust=0.5)+scale_colour_manual(values=cols) + xlim(0,0.02) +  ggtitle("sfs short/ long term infections")
```

```
## Warning: Removed 2 rows containing non-finite values (stat_density).
```

```
## Warning: Removed 6 rows containing non-finite values (stat_density).
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-1.png) 

```r
ggplot(subset(afcfs_m,value>0.01),aes(x=value,group=Var2,colour=Var2)) + geom_density(adjust=0.8) + xlim(0,0.1)+scale_colour_manual(values=cols) + ggtitle("sfs  short/ long term infections (af > 1pc)")
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-2.png) 



```r
long_term <- read.table("qstages.C1.5m.L0.01mb.G64.mtp.allgen.txt",header=T,sep="\t")
short_term <- read.table("qstages.C100.0m.L0.01mb.G64.mtp.allgen.txt",header=T,sep="\t")

lt_2pc <- cbind(long_term[,1:3],apply(long_term[,4:23],1,FUN=sum))
colnames(lt_2pc) <- c("gen","N","pc1","pc2")
lt_2pc$set="long_term"
st_2pc <- cbind(short_term[,1:3],apply(short_term[,4:23],1,FUN=sum))
colnames(st_2pc) <- c("gen","N","pc1","pc2")
st_2pc$set="short_term"

detectable <- melt(rbind(st_2pc,lt_2pc),id.vars = c("gen","set"))

cols<-c("long_term"="blue","short_term"="red")

ggplot(detectable,aes(x=gen,y=value,colour=set,group=set)) + geom_point() + facet_grid(variable ~ .,scale="free_y")  + scale_y_log10() +scale_colour_manual(values=cols)
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2-1.png) 

```r
st_melt <- melt(short_term,id.vars = c("gen","size"))
st_melt$variable <- as.numeric(gsub("X","",as.character(st_melt$variable)))
st_melt$set="short_term"
lt_melt <- melt(long_term,id.vars = c("gen","size"))
lt_melt$variable <- as.numeric(gsub("X","",as.character(lt_melt$variable)))
lt_melt$set="long_term"

cf <- rbind(subset(st_melt,gen==14),subset(lt_melt,gen==98))
```
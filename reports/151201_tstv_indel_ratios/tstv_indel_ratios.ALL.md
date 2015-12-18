library(ggplot2)
library(reshape2)
library(knitr)

```r
opts_chunk$set(fig.width=10, fig.height=9,fig.path="figure_ALL/")
```


```r
reorder <- function(M,new_order) {
  M[lower.tri(M)] = t(M)[lower.tri(M)]
  M <- M[new_order,new_order]
#  M[lower.tri(M)] <- NA
  M
}

sym <- function(M) {
  M[lower.tri(M)] = t(M)[lower.tri(M)]
  M
}


ggcolour <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}


tree_order <- c(
"Th166.12", "Th246.13", "Th245.13", "Th211.13" ,"Th092.13" ,
"Th086.07", "Th106.09",  "Th230.12","Th074.13", "Th132.11","Th162.12","Th196.12", "Th106.11", "Th117.11", "Th134.11",
"Th068.12", "Th061.13", "Th095.13"
)
clades <- c(
  rep(2,5),
  rep(3,10),
  rep(1,3))
names(clades)<-tree_order
clades
```

```
## Th166.12 Th246.13 Th245.13 Th211.13 Th092.13 Th086.07 Th106.09 Th230.12 
##        2        2        2        2        2        3        3        3 
## Th074.13 Th132.11 Th162.12 Th196.12 Th106.11 Th117.11 Th134.11 Th068.12 
##        3        3        3        3        3        3        3        1 
## Th061.13 Th095.13 
##        1        1
```


```r
tstv <- read.table("Thies_all_manual.PASS.Cls.miss0.5.LMRG.HAP.MKSNGL.tstv.tab.txt",header=T,row.names=1,sep="\t")
tstv <- reorder(tstv,tree_order)
tstv[lower.tri(tstv)] <- NA
diag(tstv)<-0

indel <- read.table("Thies_all_manual.PASS.Cls.miss0.5.LMRG.HAP.MKSNGL.in-del.tab.txt",header=T,row.names=1,sep="\t")
indel <- reorder(indel,tree_order)
diag(indel) <- 0

indelsnp <- read.table("Thies_all_manual.PASS.Cls.miss0.5.LMRG.HAP.MKSNGL.indel-snp.tab.txt",header=T,row.names=1,sep="\t")
indelsnp <- reorder(indelsnp,tree_order)
indelsnp[lower.tri(indelsnp)] <- NA

taat <- read.table("Thies_all_manual.PASS.Cls.miss0.5.LMRG.HAP.MKSNGL.taat.tab.txt",header=T,row.names=1,sep="\t")
```

```
## Warning in file(file, "rt"): cannot open file
## 'Thies_all_manual.PASS.Cls.miss0.5.LMRG.HAP.MKSNGL.taat.tab.txt': No such
## file or directory
```

```
## Error in file(file, "rt"): cannot open the connection
```

```r
taat <- reorder(taat,tree_order)
taat[lower.tri(taat)] <- NA

total_dist <- read.table("Thies_all_manual.PASS.Cls.miss0.5.LMRG.HAP.MKSNGL.dist.tab.txt",header=T,row.names=1,sep="\t")
total_dist <- reorder(total_dist,tree_order)
total_dist[lower.tri(total_dist)] <- NA
```




```r
vartypes <- merge(
  merge(
      merge(   
        melt(as.matrix(taat),value.name ="taat",factorsAsStrings = T),
        melt(as.matrix(tstv),value.name ="tstv",factorsAsStrings = T),
        by=c("Var1","Var2")
      ),
      merge(   
        melt(as.matrix(indel),value.name ="indel",factorsAsStrings = T),
        melt(as.matrix(indelsnp),value.name ="indelsnp",factorsAsStrings = T),
        by=c("Var1","Var2")
      ),
      by=c("Var1","Var2")
    ),
    melt(as.matrix(total_dist),value.name ="total",factorsAsStrings = T),
    by=c("Var1","Var2")
    )

head(vartypes)
```

```
##       Var1     Var2 taat tstv indel indelsnp total
## 1 Th061.13 Th061.13   NA 0.00  0.00       NA     0
## 2 Th061.13 Th068.12   NA   NA  1.17       NA    NA
## 3 Th061.13 Th074.13   NA   NA  1.08       NA    NA
## 4 Th061.13 Th086.07   NA   NA  1.09       NA    NA
## 5 Th061.13 Th092.13   NA   NA  1.03       NA    NA
## 6 Th061.13 Th095.13 0.54 2.92  1.00     2.44   282
```

```r
#vartypes <- vartypes[vartypes$Var1!=vartypes$Var2,]
#vartypes <- vartypes[!is.na(vartypes$total),]

vartypes$clade1 <- clades[vartypes$Var1]
vartypes$clade2 <- clades[vartypes$Var2]
vartypes$clades <- paste(clades[vartypes$Var1],clades[vartypes$Var2])

vartypes$interval_yrs <- abs(as.numeric(substr(as.character(vartypes$Var1),7,8)) - as.numeric(substr(as.character(vartypes$Var2),7,8)))
vartypes$related[vartypes$clade1 != vartypes$clade2] <- "UNRELATED"
vartypes$related[vartypes$clade1 == vartypes$clade2] <- "RELATED"

#vartypes$IBD <- vartypes$total < 2100

vartypes$Var1 <- factor(vartypes$Var1,levels=rev(tree_order))
vartypes$Var2 <- factor(vartypes$Var2,levels=rev(tree_order))

head(vartypes)
```

```
##       Var1     Var2 taat tstv indel indelsnp total clade1 clade2 clades
## 1 Th061.13 Th061.13   NA 0.00  0.00       NA     0      1      1    1 1
## 2 Th061.13 Th068.12   NA   NA  1.17       NA    NA      1      1    1 1
## 3 Th061.13 Th074.13   NA   NA  1.08       NA    NA      1      3    1 3
## 4 Th061.13 Th086.07   NA   NA  1.09       NA    NA      1      3    1 3
## 5 Th061.13 Th092.13   NA   NA  1.03       NA    NA      1      2    1 2
## 6 Th061.13 Th095.13 0.54 2.92  1.00     2.44   282      1      1    1 1
##   interval_yrs   related
## 1            0   RELATED
## 2            1   RELATED
## 3            0 UNRELATED
## 4            6 UNRELATED
## 5            0 UNRELATED
## 6            0   RELATED
```

all matrices

```r
relcol <- scale_color_manual(values=c("white","black"))
vxlab <- theme(axis.text.x = element_text(angle = 90, hjust = 1))  

#total pairwise distance (DISCORDS)
ggplot(subset(vartypes,!is.na(total)),aes(x=Var1,y=Var2,fill=total,label=total,colour=related)) + geom_tile() + scale_fill_gradient(trans="log") + geom_text(size=4) + relcol + vxlab
```

![plot of chunk unnamed-chunk-4](figure_ALL/unnamed-chunk-4-1.png) 

```r
#ggplot(subset(vartypes,!is.na(total)),aes(x=Var1,y=Var2,fill=total,label=total,colour=related)) + geom_tile() + scale_fill_gradient(trans="log") + geom_text(size=4) + relcol + vxlab

#Ts:Tv ratio
ggplot(subset(vartypes,!is.na(tstv)),aes(x=Var1,y=Var2,fill=tstv,label=tstv)) + geom_tile() + scale_fill_gradient(trans="log") + geom_text(size=4,colour="white") + relcol + vxlab
```

![plot of chunk unnamed-chunk-4](figure_ALL/unnamed-chunk-4-2.png) 

ts:tv / in:del / distance scatterplots 

```r
#Ts:Tv x dist
ggplot(subset(vartypes,tstv > 0),aes(x=total,y=tstv,colour=related)) + geom_point() + vxlab + ylim(0,5)
```

![plot of chunk unnamed-chunk-5](figure_ALL/unnamed-chunk-5-1.png) 

```r
#ggplot(subset(vartypes,tstv > 0),aes(x=total,y=tstv,colour=related)) + geom_point() + vxlab + scale_y_log10(limits=c(0.2,5))
#IN:DEL x dist
ggplot(subset(vartypes,tstv > 0),aes(x=total,y=indel,colour=related)) + geom_point() + vxlab + scale_y_log10(limits=c(0.5,2))
```

![plot of chunk unnamed-chunk-5](figure_ALL/unnamed-chunk-5-2.png) 

```r
#Ts:Tv x IN:DEL
#ggplot(subset(vartypes,tstv > 0),aes(x=indel,y=tstv,colour=related,size=1/total)) + geom_point() + vxlab
ggplot(subset(vartypes,tstv > 0),aes(x=indel,y=tstv,colour=related,size=1/total)) + geom_point() + vxlab + scale_x_log10(limits=c(0.5,2))
```

![plot of chunk unnamed-chunk-5](figure_ALL/unnamed-chunk-5-3.png) 

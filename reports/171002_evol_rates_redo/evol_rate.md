

```r
library(ape)
library(adegenet)
library(knitr)
library(igraph)
library(RColorBrewer)
library(gridExtra)
library(ggplot2)
opts_chunk$set(fig.width=9, fig.height=9)
opts_chunk$set(dev=c('png','postscript'))
```



```r
sym <- function(M) {
  M[lower.tri(M)] = t(M)[lower.tri(M)]
  M
}

makeNet <- function(distance_matrix, meta_file, ngroups=3) {
  if (class(distance_matrix)=="dist") {
    D = distance_matrix
    mat <- as.matrix(D)
  } else {
    mat <- read.table(distance_matrix,sep="\t")
    D <- as.dist(sym(mat))
    }
  
  clust <- gengraph(D,ngrp=3)
  names <- colnames(mat)
  mat <- as.matrix(mat)

  meta <- read.table("Thies_metadata_1701.txt",sep="\t",header=T)
  colnames(meta)[1]<-"name"
  meta <- meta[!is.na(meta$Age),]

  rownames(meta) <- meta$name
  meta <- meta[names,]
  coll <- as.Date(as.character(meta$Date),"%d/%m/%Y",origin = "2000-01-01")
  #coll <- as.Date(paste("1","jan",meta$year,sep=""),"%d%b%Y")
  names(coll)<-meta$name
  meta$year <- as.numeric(format(coll,'%Y'))+2000

  
#   meta <- read.table("daniels.thies.CA.txt",sep="\t",header=T)
#   rownames(meta) <- meta$name
#   meta <- meta[names,]
#   coll <- as.Date(paste("1","jan",meta$year,sep=""),"%d%b%Y")
#   names(coll)<-meta$name

  name1 <- names[clust$clust$membership==1]
  name2 <- names[clust$clust$membership==2]
  name3 <- names[clust$clust$membership==3]
  year1 <- meta$year[clust$clust$membership==1]
  year2 <- meta$year[clust$clust$membership==2]
  year3 <- meta$year[clust$clust$membership==3]
  coll1 <- coll[name1]
  coll2 <- coll[name2]
  coll3 <- coll[name3]
  dist1 <- mat[name1,name1]
  dist2 <- mat[name2,name2]
  dist3 <- mat[name3,name3]
  
  res1 <- seqTrack(dist1, x.names=name1, x.dates=coll1)
  res2 <- seqTrack(dist2, x.names=name2, x.dates=coll2)
  res3 <- seqTrack(dist3, x.names=name3, x.dates=coll3)

  res1$year <- year1
  res2$year <- year2
  res3$year <- year3
  res1$name <- name1
  res2$name <- name2
  res3$name <- name3
  
  list(res1,res2,res3)
}
```

```r
ifile <- "Thies_all_manual.PASS.Cls.miss0.5.LMRG.HAP.INDEL.recode.vcf.dist.tab.txt"
sfile <- "Thies_all_manual.PASS.Cls.miss0.5.LMRG.HAP.SNP.recode.vcf.dist.tab.txt"
corefile <- "Thies_all_manual.PASS.Cls.miss0.5.LMRG.HAP.CALLBOTH.vcf.dist.tab.txt"
distmatS <- as.dist(sym(read.table(sfile,sep="\t")))
distmatI <- as.dist(sym(read.table(ifile,sep="\t")))
distmatIS <- distmatS+distmatI
distmatCore <- as.dist(sym(read.table(corefile,sep="\t")))

indelNets <- makeNet(distmatI)
snpNets <- makeNet(distmatS)

varNets <- makeNet(distmatIS)
```




```r
printGraph <- function(graph,colours,title) {
  cols <- brewer.pal((max(graph$year)-min(graph$year))+1, colours)

  ts=1 #textsize
  ig <- as.igraph(graph)
  tree <- layout_as_tree(ig,flip.y = F)[,c(2,1)]
  V(ig)$name <- graph$name
  V(ig)$color <- cols[graph$year-min(graph$year)+1]
  V(ig)$label.cex <- ts
 #frame()
  plot(ig,layout=tree,main=title,vertex.size=25,
     edge.color="black",edge.label.cex=1.5,edge.label.family="Helvetica",vertex.label.family="Helvetica")

}

#printGraph(indelNets[[1]],"Greens","clade 2, indel graph")

#printGraph(snpNets[[1]],"Greens","clade 2, indel graph")

#printGraph(varNets[[2]],"Blues","clade 29 all vars")
printGraph(varNets[[1]],"Greens","clade 29 all vars")
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-1.png)

```r
#printGraph(varNets[[3]],"Reds","clade 29 all vars")
```




```r
net <- varNets[[1]]

indDists <- as.matrix(distmatI)
snpDists <- as.matrix(distmatS)
coreDists <- as.matrix(distmatCore)
net$direct<-T
netAll <- net
getAnces <- function(net, leaf, ids=character()) {
    id = ids[1]
    ances = netAll[id,"ances"]
    date = netAll[id,"date"]
    ances.date = netAll[id,"ances.date"]
#    if(is.na(steps)) {steps = netAll[id,"steps"]} else {write(steps,stderr())}
    #write(paste(leaf,id,ances,sep="\t"),stderr())
    if (!is.na(ances)) {
      if (id != leaf) {
        leaf.date = net[leaf,"date"]
        leaf.year = net[leaf,"year"]
        leaf.name = net[leaf,"name"]
        nextI = dim(net)[[1]]+1
#        write(paste("  adding",nextI,":",leaf,ances,sep=" "),stderr())
        net[nextI,1:2] <- c(leaf, ances)
        net[nextI,4] <- leaf.date
        net[nextI,5] <- ances.date
      
        net[nextI,6] <- leaf.year
        net[nextI,7] <- leaf.name
        net[nextI,8] <- F
        }
      net <- getAnces(net,leaf,c(ances))
    }
    if (length(ids) > 1){
       ids <- ids[ids != id]
#       write(paste("proceeding:",ids[1],sep=" "),stderr())
       net <- getAnces(net,ids[1],ids)
#       steps <- steps+1
#       net <- getAnces(net,ids[1],ids,steps)
  }
net
}
getAnces(net,net$id[1],net$id)
```

```
##          id ances weight       date ances.date year     name direct
## Th086.07  1    NA     NA 0007-10-11       <NA> 2007 Th086.07   TRUE
## Th106.09  2     1    576 0009-10-26 0007-10-11 2009 Th106.09   TRUE
## Th106.11  3     2    358 0011-10-27 0009-10-26 2011 Th106.11   TRUE
## Th117.11  4     3    111 0011-11-02 0011-10-27 2011 Th117.11   TRUE
## Th132.11  5     2    121 0011-11-11 0009-10-26 2011 Th132.11   TRUE
## Th134.11  6     3    105 0011-11-14 0011-10-27 2011 Th134.11   TRUE
## Th162.12  7     2    125 0012-10-10 0009-10-26 2012 Th162.12   TRUE
## Th196.12  8     5    139 0012-10-25 0011-11-11 2012 Th196.12   TRUE
## Th230.12  9     5    106 0012-11-09 0011-11-11 2012 Th230.12   TRUE
## Th074.13 10     2     88 0013-09-23 0009-10-26 2013 Th074.13   TRUE
## 11        3     1     NA 0011-10-27 0007-10-11 2011 Th106.11  FALSE
## 12        4     2     NA 0011-11-02 0009-10-26 2011 Th117.11  FALSE
## 13        4     1     NA 0011-11-02 0007-10-11 2011 Th117.11  FALSE
## 14        5     1     NA 0011-11-11 0007-10-11 2011 Th132.11  FALSE
## 15        6     2     NA 0011-11-14 0009-10-26 2011 Th134.11  FALSE
## 16        6     1     NA 0011-11-14 0007-10-11 2011 Th134.11  FALSE
## 17        7     1     NA 0012-10-10 0007-10-11 2012 Th162.12  FALSE
## 18        8     2     NA 0012-10-25 0009-10-26 2012 Th196.12  FALSE
## 19        8     1     NA 0012-10-25 0007-10-11 2012 Th196.12  FALSE
## 20        9     2     NA 0012-11-09 0009-10-26 2012 Th230.12  FALSE
## 21        9     1     NA 0012-11-09 0007-10-11 2012 Th230.12  FALSE
## 22       10     1     NA 0013-09-23 0007-10-11 2013 Th074.13  FALSE
```

```r
getEvolRates <- function(net,mat,ng,ances=TRUE) {
  if(ances) {netAll <- getAnces(net, net$id[1],net$id)}
  else {netAll <- net}
  
  D <- as.dist(sym(mat))
  names <- colnames(mat)
  mat <- as.matrix(sym(mat))
  outgroups <- setdiff(names,net$name)
  
  outtab <- data.frame(from=character(),
                 to=character(),
                 out=character(),
                 time=integer(),
                 years=integer(),
                 distance=numeric(),
                 rate=numeric(),
             #    direct=logical(),
                 stringsAsFactors=FALSE)
  
  for (i in c(1:dim(netAll)[[1]])) {
    if(!is.na(netAll$ances[i])) {
      sample = netAll$name[i]
      ances = netAll$name[netAll$ances[i]]
      #isdirect <- netAll$direct[i]
      dists = c()
      
      #time elapsed in years
      #time <- as.integer(netAll$year[i] - netAll$year[netAll$ances[i]])
      years <- as.integer(netAll$year[i] - netAll$year[netAll$ances[i]])
      
      #time elapsed in months
      time <- round(as.integer(netAll$date[i]-netAll$ances.date[i])*12/365,2)
      
      for (o in outgroups) {
        dist <- mat[sample,o]-mat[ances,o]
        dists <- c(dists, dist)
        outtab[dim(outtab)[1]+1,1:3]=c(ances,sample,o)
        outtab[dim(outtab)[1],4:7]=c(time,years,dist,round(dist/time,2))  
        #outtab[dim(outtab)[1],8]=isdirect
      }
      dists <- round(dists/time)
      #paste(dists,collapse = ","),
      #write(paste(ances,sample,netAll$date[i],netAll$ances.date[i],time,round(mean(dists),2),round(sd(dists),2),sep="\t"),stdout())
#      write(paste(ances,sample,outgroups,time,dists,sep="\t"),stdout())
    }
  }
  outtab$pair <- paste(outtab$from,outtab$to)
  outtab
}
test <- getEvolRates(varNets[[1]],indDists)


g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
```


```r
opts_chunk$set(fig.width=18, fig.height=8)
opts_chunk$set(dev=c('png','postscript'))
```






```r
ir1 <- getEvolRates(varNets[[1]],indDists,ances=F)
sr1 <- getEvolRates(varNets[[1]],snpDists,ances=F)  
isr1 <- getEvolRates(varNets[[1]],indDists+snpDists,ances=F) 
ir1 <- subset(ir1,years>0)
sr1 <- subset(sr1,years>0)
isr1 <- subset(isr1,years>0)

#INDELs ONLY
iplot <- ggplot(subset(ir1,years>0),aes(x=years,y=rate)) + 
  geom_boxplot(aes(group=years),fill=NA) + geom_hline(aes(yintercept=mean(rate,na.rm=T)),linetype=3) + 
  geom_point(aes(colour=to,shape=from),position=position_jitter(width=0.5),size=3)+
  geom_label(aes(label=paste(round(mean(ir1$rate,na.rm=T),2)," (+/- ",round(sd(ir1$rate,na.rm=T),1),")",sep=""),
                 y=mean(ir1$rate,na.rm=T), x=max(years)),size=5)+
  ggtitle("indel rate, clade 29") + ylab("INDELs/mth/gen") + xlab("years")  + ylim(-2,9) + #xlim(0,7.5) +
  theme(legend.position="none")
#SNPs ONLY
splot <- ggplot(subset(sr1,years>0),aes(x=years,y=rate)) + 
  geom_boxplot(aes(group=years),fill=NA) + geom_hline(aes(yintercept=mean(rate,na.rm=T)),linetype=3) + 
  geom_point(aes(colour=to,shape=from),position=position_jitter(width=0.5),size=3)+
  geom_label(aes(label=paste(round(mean(sr1$rate,na.rm=T),2)," (+/- ",round(sd(sr1$rate,na.rm=T),1),")",sep=""),
                 y=mean(sr1$rate,na.rm=T)+0.8, x=max(years)),size=5)+
  ggtitle("substitution rate, clade 29") + ylab("SNPs/mth/gen") + xlab("years") + ylim(-2,9) #+ xlim(0,7.5)
#splot
legend <- g_legend(splot)
grid.arrange(iplot, splot + theme(legend.position="none"),legend, nrow=1,widths=c(10,10,1),heights=c(1))
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-1.png)


```r
corer1 <- getEvolRates(varNets[[1]],coreDists,ances=F) 
inaccr1 <- getEvolRates(varNets[[1]],(indDists+snpDists)-coreDists,ances=F) 
corer1 <- subset(corer1,years>0)
inaccr1 <- subset(inaccr1,years>0)


#CORE ONLY
coreplot <- ggplot(subset(corer1,years>0),aes(x=years,y=rate)) + 
  geom_boxplot(aes(group=years),fill=NA) + geom_hline(aes(yintercept=mean(rate,na.rm=T)),linetype=3) + 
  geom_point(aes(colour=to,shape=from),position=position_jitter(width=0.5),size=3)+
  geom_label(aes(label=paste(round(mean(corer1$rate,na.rm=T),2)," (+/- ",round(sd(corer1$rate,na.rm=T),1),")",sep=""),
                 y=mean(corer1$rate,na.rm=T)+0.2,x=max(years)),size=5)+
  ggtitle("evol rate, core") + ylab("INDELs/mth/gen") + xlab("years")  + ylim(-2,9) + #xlim(0,7.5) +
  theme(legend.position="none")
#INACCESSIBLE ONLY
inaccplot <- ggplot(subset(inaccr1,years>0),aes(x=years,y=rate)) + 
  geom_boxplot(aes(group=years),fill=NA) + geom_hline(aes(yintercept=mean(rate,na.rm=T)),linetype=3) + 
  geom_point(aes(colour=to,shape=from),position=position_jitter(width=0.5),size=3)+
  geom_label(aes(label=paste(round(mean(inaccr1$rate,na.rm=T),2)," (+/- ",round(sd(inaccr1$rate,na.rm=T),1),")",sep=""),
                 y=mean(inaccr1$rate,na.rm=T)+1.0,x=max(years)),size=5)+
  ggtitle("evol rate, inaccessible") + ylab("SNPs/mth/gen") + xlab("years") + ylim(-2,9) #+ xlim(0,7.5)
#splot
legend <- g_legend(coreplot)
```

```
## Warning: Removed 2 rows containing non-finite values (stat_boxplot).
```

```
## Warning: Removed 2 rows containing missing values (geom_point).
```

```
## Error in tmp$grobs[[leg]]: attempt to select less than one element
```

```r
grid.arrange(coreplot, inaccplot + theme(legend.position="none"),legend, nrow=1,widths=c(10,10,1),heights=c(1))
```

```
## Warning: Removed 2 rows containing non-finite values (stat_boxplot).

## Warning: Removed 2 rows containing missing values (geom_point).
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-1.png)



```r
# 
# #INDELs ONLY
# iplot <- ggplot(subset(ir1,years>0),aes(x=years,y=rate)) + 
#   geom_boxplot(aes(group=years),fill=NA) + geom_hline(aes(yintercept=mean(rate,na.rm=T)),linetype=3) + 
#   geom_point(aes(colour=direct,shape=to),position=position_jitter(width=0.5),size=3)+
#   geom_label(aes(label=paste(round(mean(rate,na.rm=T),2)," (+/- ",round(sd(rate,na.rm=T),1),")",sep=""),y=mean(rate,na.rm=T),x=7.3),size=5)+
#   ggtitle("indel rate, clade 2") + ylab("INDELs/mth/gen") + xlab("years")  + ylim(-40,9) + #xlim(0,7.5) +
#   theme(legend.position="none")
# #SNPs ONLY
# splot <- ggplot(subset(sr1,years>0),aes(x=years,y=rate)) + 
#   geom_boxplot(aes(group=years),fill=NA) + geom_hline(aes(yintercept=mean(rate,na.rm=T)),linetype=3) + 
#   geom_point(aes(colour=direct,shape=to),position=position_jitter(width=0.5),size=3)+
#   geom_label(aes(label=paste(round(mean(rate,na.rm=T),2)," (+/- ",round(sd(rate,na.rm=T),1),")",sep=""),y=mean(rate,na.rm=T),x=max(years)),size=5)+
#   ggtitle("substitution rate, clade 2") + ylab("SNPs/mth/gen") + xlab("years") + ylim(-40,9) #+ xlim(0,7.5)
# #splot
# legend <- g_legend(splot)
# grid.arrange(iplot, splot + theme(legend.position="none"),legend, nrow=1,widths=c(10,10,1),heights=c(1))
```



```r
opts_chunk$set(fig.width=12, fig.height=8)
opts_chunk$set(dev=c('png','postscript'))
```

```r
#SNPs + INDELs
ggplot(subset(isr1,years>0),aes(x=years,y=rate)) + 
  geom_boxplot(aes(group=years),fill=NA) + geom_hline(aes(yintercept=median(rate,na.rm=T)),linetype=3) + 
  geom_point(aes(colour=to,shape=from),position=position_jitter(width=0.5),size=3) +
  geom_label(aes(label=paste(round(mean(isr1$rate,na.rm=T),2)," (+/- ",round(sd(isr1$rate,na.rm=T),1),")",sep=""),
                 y=mean(isr1$rate,na.rm=T)+0.1,x=max(years)),size=5)+
  ggtitle("evolutionary rate, clade 29") + ylab("mut/mth/gen") + xlab("years")  #+ ylim(-2,9) + xlim(0,7.5)
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-11-1.png)



```r
#get mean rates for > 0 year, in years
irate <- round(mean(ir1$rate)*12,3)
srate <- round(mean(sr1$rate)*12,3)
israte <- round(mean(isr1$rate)*12,3)
corerate <- round(mean(corer1$rate)*12,3)
inaccrate <- round(mean(inaccr1$rate)*12,3)
```


```r
sppsubs <- read.table("subrate_relaxed_otherspp2.txt",sep="\t",header=T,stringsAsFactors = F)
colnames(sppsubs)[4:6] <- c("gsize","subrate","subrate_genome")
sppsubs$Group <- factor(sppsubs$Group,levels=c("RNA virus","dsDNA virus","ssDNA virus","Bacteria","Eukaryote"),ordered = T)

pfmingsize <- 15180000
pfgsize<-18446000
pfiagsize <- 18446000-15180000

# sppsubs <- rbind(sppsubs,
#       c("P.falciparum","eukaryote",pfgsize,signif(israte/pfgsize,2),israte,"among",8,"ours")
#       )

sppsubs$vartype<-"SNP"
sppsubs$region<-"all"
sppsubs <- rbind(sppsubs,
      c("P.falciparum","P.falciparum (snp)","Eukaryote",pfgsize,signif(srate/pfgsize,2),srate,"among",8,"ours","SNP","all"),
      c("P.falciparum","P.falciparum (ind)","Eukaryote",pfgsize,signif(irate/pfgsize,2),irate,"among",8,"ours","INDEL","all"),
      c("P.falciparum","P.falciparum (all)","Eukaryote",pfgsize,signif(israte/pfgsize,2),israte,"among",8,"ours","SNP+INDEL","all"),
      c("P.falciparum","P.falciparum (core)","Eukaryote",pfmingsize,signif(corerate/pfmingsize,2),corerate,"among",8,"ours","SNP+INDEL","accessible"),
      c("P.falciparum","P.falciparum (inacc)","Eukaryote",pfgsize - pfmingsize,signif(inaccrate/pfiagsize,2),inaccrate,"among",8,"ours","SNP+INDEL","inaccessible")
      )

sppsubs$subrate <- as.numeric(sppsubs$subrate)
sppsubs$subrate_genome <- as.numeric(sppsubs$subrate_genome)
sppsubs$gsize <- as.numeric(sppsubs$gsize)
sppsubs$vartype <- factor(sppsubs$vartype,levels=c("SNP","INDEL","SNP+INDEL"),ordered=T)
#sppsubs$Label <- factor(sppsubs$Label,levels=sppsubs$Label[order(sppsubs$Group,sppsubs$subrate_genome)],ordered=T)


pathOrder <- c("HCV","EBOV", "FMDV",   "HIV-M","HFLUV", "RABV",  
 "H. pylori", "S. aureus", "S. pneumoniae", "V. cholerae", "M. tuberculosis", 
"P.falciparum")

sppsubs$Pathogen <- factor(sppsubs$Pathogen,levels=pathOrder,ordered=T)

gsizes <- c(1e3,1e4,1e5,1e6,1e7,1e8)
rates <- data.frame("gsize"=gsizes,"years"=1/gsizes,"months"=12/gsizes,"weeks"=52/gsizes)
ratesM <- melt(rates,id.vars = "gsize",value.name ="subrate")
#ratesM <- ratesM[ratesM$gsize>0,]
floors <- data.frame("gsize"=c(rep(max(ratesM$gsize),3),rep(min(ratesM$gsize),3)),
      "variable"=unique(as.character(ratesM$variable)),
      "subrate"=rep(0,6))

ratesM <- rbind(ratesM,floors)
ratesM$variable <- factor(ratesM$variable,levels=c("weeks","months","years"),ordered=T)

fcol <- c("#AAAAAA","#CCCCCC","#EEEEEE")
names(fcol) = c("years","months","weeks")

dotsubs <- c("H. pylori", "M. tuberculosis", "S. aureus (MRSA)", "S. pneumoniae", 
"V. cholerae", "HCV (Hep C)", "EBOV (Ebola)", "HIV-M", # "HFLUV (H1N1 2009)", 
"HFLUV (H1N1)", "RABV", "FMDV",    "P.falciparum (snp)", # "P.falciparum (ind)", 
"P.falciparum (all)", "P.falciparum (core)")

ggplot(subset(sppsubs, Label %in% dotsubs),aes(x=gsize,y=subrate),) +  
  geom_polygon(data=ratesM,aes(group=variable,fill=variable)) +
  scale_fill_manual(values = fcol) + 
  geom_label(aes(label=Label,colour=Group),nudge_x = -0.2,hjust=1) + 
  geom_point(aes(label=Pathogen,colour=Group,shape=vartype),size=5) + 
  scale_y_log10() + scale_x_log10(limits=c(1e3,1e8),breaks=c(1e3,1e4,1e5,1e6,1e7,1e8)) +
  xlab("Genome Size (bp)")+ylab("mutation rate (muts / site / year)")+
  theme(panel.background = element_rect(fill = F))
```

```
## Warning: Ignoring unknown aesthetics: label
```

```
## Warning: Transformation introduced infinite values in continuous y-axis
```

![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-13-1.png)

```r
#dput(sppsubs$Label)

yearRate <- 1/(sppsubs$subrate*sppsubs$gsize)
monthRate <- 1/(sppsubs$subrate*sppsubs$gsize/12)
weekRate <- 1/(sppsubs$subrate*sppsubs$gsize/52)

period <- rep("yrs",length(yearRate))
newmut <- round(yearRate,1)
period[yearRate < 1] <- "mths"
newmut[yearRate < 1] <- round(monthRate,1)[yearRate < 1]
period[monthRate < 1] <- "wks"
newmut[monthRate < 1] <- round(weekRate,1)[monthRate < 1]
sppsubs$newmut <- paste(newmut,period,sep=" ")

barSamps <- c("M. tuberculosis", "S. aureus (MRSA)","H. pylori",  
"EBOV (Ebola)", "HIV-M", 
"HFLUV (H1N1)", "P.falciparum (inacc)",
"P.falciparum (core)")
labelSamps <- c("M. tuberculosis", "S. aureus (MRSA)","H. pylori", 
"EBOV (Ebola)", "HIV-M", 
"HFLUV (H1N1)", "P.falciparum (all)")
midLabelSamps <- c("P.falciparum (inacc)",
"P.falciparum (core)")

# ggplot(subset(sppsubs,Label %in% barSamps),aes(x=Pathogen,y=subrate * gsize),) +  
#   geom_bar(aes(fill=Group,alpha=region),stat = "identity",size=3) +
#   geom_label(data=subset(sppsubs,Label %in% labelSamps),aes(label=newmut))+
#   geom_label(data=subset(sppsubs,Label %in% midLabelSamps),aes(label=newmut, y=(subrate * gsize)/2))+
# #  scale_alpha_manual(values=c("SNP"=0.9,"INDEL"=0.6,"SNP+INDEL"=1))+
#   scale_alpha_manual(values=c("all"=1,"accessible"=1,"inaccessible"=0.8))+
#   ylab("mutation rate (muts / genome / year)")+
#   theme(panel.background = element_rect(fill = F),axis.text.x = element_text(angle = 90, hjust = 1))
# 

write.table(sppsubs,stdout(),sep="\t",quote=F)
```

```
## Pathogen	Label	Group	gsize	subrate	subrate_genome	Host.scale	Date.range..yrs.	Source	vartype	region	newmut
## 1	H. pylori	H. pylori	Bacteria	1600000	1.38e-05	22	among	15	[S5]	SNP	all	2.4 wks
## 2	M. tuberculosis	M. tuberculosis	Bacteria	4400000	6.8e-08	0.3	among	10	[S9]	SNP	all	3.3 yrs
## 3	S. aureus	S. aureus (MRSA)	Bacteria	2900000	3.3e-06	9.57	among	25	[S13]	SNP	all	1.3 mths
## 4	S. pneumoniae	S. pneumoniae	Bacteria	2200000	1.6e-06	3.52	among	25	[S15]	SNP	all	3.4 mths
## 5	V. cholerae	V. cholerae	Bacteria	4e+06	8.3e-07	3.32	among	60	[S16]	SNP	all	3.6 mths
## 6	HCV	HCV (Hep C)	RNA virus	9000	0.001	9	among	35	[S20]	SNP	all	1.3 mths
## 7	EBOV	EBOV (Ebola)	RNA virus	19000	0.002	9	among	35	[S20]	SNP	all	1.4 wks
## 8	HIV-M	HIV-M	RNA virus	10000	0.0033	33	among	50	[S21]	SNP	all	1.6 wks
## 9	HFLUV	HFLUV (H1N1 2009)	RNA virus	13600	0.0036	48.96	among	1	[S22]	SNP	all	1.1 wks
## 10	HFLUV	HFLUV (H1N1)	RNA virus	13600	0.002	27.2	among	60	[S23]	SNP	all	1.9 wks
## 11	RABV	RABV	RNA virus	12000	2e-04	2.4	among	22	[S24], Biek, unpubl.	SNP	all	5 mths
## 12	FMDV	FMDV	RNA virus	8100	0.0013	9	among	35	[S20]	SNP	all	1.1 mths
## 13	P.falciparum	P.falciparum (snp)	Eukaryote	18446000	1.4e-06	25.734	among	8	ours	SNP	all	2 wks
## 14	P.falciparum	P.falciparum (ind)	Eukaryote	18446000	5.1e-07	9.317	among	8	ours	INDEL	all	1.3 mths
## 15	P.falciparum	P.falciparum (all)	Eukaryote	18446000	1.9e-06	35.068	among	8	ours	SNP+INDEL	all	1.5 wks
## 16	P.falciparum	P.falciparum (core)	Eukaryote	15180000	6.6e-07	10.056	among	8	ours	SNP+INDEL	accessible	1.2 mths
## 17	P.falciparum	P.falciparum (inacc)	Eukaryote	3266000	7.7e-06	25.005	among	8	ours	SNP+INDEL	inaccessible	2.1 wks
```

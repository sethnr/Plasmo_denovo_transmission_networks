
```r
library(ape)
library(adegenet)
library(knitr)
library(igraph)
library(RColorBrewer)
library(gridExtra)
library(ggmap)

opts_chunk$set(fig.width=9, fig.height=9)
opts_chunk$set(dev=c('png','postscript'))
```



```r
sym <- function(M) {
  M[lower.tri(M)] = t(M)[lower.tri(M)]
  M
}

makeNetSingle <- function(distance_matrix, meta) {
  if (class(distance_matrix)=="dist") {
    D = distance_matrix
    mat <- as.matrix(D)
  } else {
    mat <- read.table(distance_matrix,sep="\t")
    D <- as.dist(sym(mat))
    }
  
  names <- colnames(mat)
  mat <- as.matrix(mat)

  #limit meta to only these files
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

  year <- meta$year
  
  res <- seqTrack(mat, x.names=names, x.dates=coll)

  res$year <- year
  res$name <- names
  res
}
```


```r
# ifile <- "Thies_all_manual.PASS.Cls.miss0.5.LMRG.HAP.INDEL.recode.vcf.dist.tab.txt"
# sfile <- "Thies_all_manual.PASS.Cls.miss0.5.LMRG.HAP.SNP.recode.vcf.dist.tab.txt"
# corefile <- "Thies_all_manual.PASS.Cls.miss0.5.LMRG.HAP.CALLBOTH.vcf.dist.tab.txt"
# distmatS <- as.dist(sym(read.table(sfile,sep="\t")))
# distmatI <- as.dist(sym(read.table(ifile,sep="\t")))
# distmatIS <- distmatS+distmatI
# distmatCore <- as.dist(sym(read.table(corefile,sep="\t")))

# indelNets <- makeNet(distmatI)
# snpNets <- makeNet(distmatS)

# varNets <- makeNet(distmatIS)
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
     edge.color="black",edge.label.cex=1.5,edge.label.family="Arial",vertex.label.family="Arial")

}
```




```r
opts_chunk$set(fig.width=12, fig.height=9)
opts_chunk$set(dev=c('png'))
```


```r
net <- varNets[[1]]
```

```
## Error in eval(expr, envir, enclos): object 'varNets' not found
```

```r
indDists <- as.matrix(distmatI)
```

```
## Error in as.matrix(distmatI): error in evaluating the argument 'x' in selecting a method for function 'as.matrix': Error: object 'distmatI' not found
```

```r
snpDists <- as.matrix(distmatS)
```

```
## Error in as.matrix(distmatS): error in evaluating the argument 'x' in selecting a method for function 'as.matrix': Error: object 'distmatS' not found
```

```r
coreDists <- as.matrix(distmatCore)
```

```
## Error in as.matrix(distmatCore): error in evaluating the argument 'x' in selecting a method for function 'as.matrix': Error: object 'distmatCore' not found
```

```r
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
##          id ances weight       date ances.date year     name  type boot
## Th086.07  1    NA     NA 0007-10-11       <NA> 2007 Th086.07  core    0
## Th106.09  2     1    576 0009-10-26 0007-10-11 2009 Th106.09  core    0
## Th106.11  3     2    358 0011-10-27 0009-10-26 2011 Th106.11  core    0
## Th117.11  4     3    111 0011-11-02 0011-10-27 2011 Th117.11  core    0
## Th132.11  5     2    121 0011-11-11 0009-10-26 2011 Th132.11  core    0
## Th134.11  6     3    105 0011-11-14 0011-10-27 2011 Th134.11  core    0
## Th162.12  7     2    125 0012-10-10 0009-10-26 2012 Th162.12  core    0
## Th196.12  8     5    139 0012-10-25 0011-11-11 2012 Th196.12  core    0
## Th230.12  9     5    106 0012-11-09 0011-11-11 2012 Th230.12  core    0
## Th074.13 10     2     88 0013-09-23 0009-10-26 2013 Th074.13  core    0
## 11        3     1     NA 0011-10-27 0007-10-11 2011 Th106.11 FALSE   NA
## 12        4     2     NA 0011-11-02 0009-10-26 2011 Th117.11 FALSE   NA
## 13        4     1     NA 0011-11-02 0007-10-11 2011 Th117.11 FALSE   NA
## 14        5     1     NA 0011-11-11 0007-10-11 2011 Th132.11 FALSE   NA
## 15        6     2     NA 0011-11-14 0009-10-26 2011 Th134.11 FALSE   NA
## 16        6     1     NA 0011-11-14 0007-10-11 2011 Th134.11 FALSE   NA
## 17        7     1     NA 0012-10-10 0007-10-11 2012 Th162.12 FALSE   NA
## 18        8     2     NA 0012-10-25 0009-10-26 2012 Th196.12 FALSE   NA
## 19        8     1     NA 0012-10-25 0007-10-11 2012 Th196.12 FALSE   NA
## 20        9     2     NA 0012-11-09 0009-10-26 2012 Th230.12 FALSE   NA
## 21        9     1     NA 0012-11-09 0007-10-11 2012 Th230.12 FALSE   NA
## 22       10     1     NA 0013-09-23 0007-10-11 2013 Th074.13 FALSE   NA
##          direct
## Th086.07   TRUE
## Th106.09   TRUE
## Th106.11   TRUE
## Th117.11   TRUE
## Th132.11   TRUE
## Th134.11   TRUE
## Th162.12   TRUE
## Th196.12   TRUE
## Th230.12   TRUE
## Th074.13   TRUE
## 11           NA
## 12           NA
## 13           NA
## 14           NA
## 15           NA
## 16           NA
## 17           NA
## 18           NA
## 19           NA
## 20           NA
## 21           NA
## 22           NA
```

```r
getEvolRates <- function(net,mat,ng) {
  netAll <- getAnces(net, net$id[1],net$id)
  
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
```

```
## Error in getAnces(net, net$id[1], net$id): object 'varNets' not found
```

```r
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
ped <- "Thies_all_manual.PASS.Cls.miss0.5.LMRG.HAP.DISCORDS.vcf.gz_NJnex/Thies_all_manual.PASS.Cls.miss0.5.LMRG.HAP.DISCORDS.gz.ped"
genos <- read.table(ped,colClasses="character")
inds <- genos[,1]
genos <- genos[,seq(7,dim(genos)[[2]],2)]
rownames(genos)=inds


  meta <- read.table("Thies_metadata_1701.txt",sep="\t",header=T)
  colnames(meta)[1]<-"name"
  meta <- meta[!is.na(meta$Age),]


  #tab<-"./Thies_all_manual.PASS.Cls.miss0.5.LMRG.HAP.alleles.tab"
  fileroot<-"ThiesDiscoAlleles_ped"
  
  #read in ped file and get hamming distances:
  #alleleTab <- read.table(tab,colClasses="character",header=T,na.strings = c("."))
  
  c1 <- c("Th086.07","Th106.09","Th106.11","Th117.11","Th132.11","Th134.11","Th162.12","Th196.12","Th230.12","Th074.13")
  #genos <- t(data.matrix(alleleTab[3:dim(alleleTab)[2]]))


  genos <- genos[c1,]
  inds <- c1


#  inds<-row.names(genos)
#   DNdistmat = matrix(nrow=length(c1),ncol=length(c1))
#   colnames(DNdistmat) = c1
#   rownames(DNdistmat) = c1
# 
#   outs <- inds[(!inds %in% c1)]
#   write("calculating de novo distance matrix",stderr())
# 
#   for (i in c1){
#       for (j in c1){
#         for (o in outs){
#             if(length(filled) > 0) {
#             meta[meta$name==i,"Date"] > meta[meta$name==j,"Date"]
#             distmat[i,j] = 
#         }else {
#           distmat[i,j]=-1
#         }
#       }
#   }
# 
  


  distmat = matrix(nrow=length(inds),ncol=length(inds))
  colnames(distmat) = inds
  rownames(distmat) = inds

  write("calculating distance matrix",stderr())
  for (i in inds){
      for (j in inds){
      filled = intersect(which(genos[i,] != 0), which(genos[j,] !=0))
      write(paste("calculating",i,"v",j),stderr())
      write(length(filled),stderr())
  
  #    distmat[i,j] = sum(genos[i,]!=genos[j,])
      if(length(filled) > 0) {
          distmat[i,j] = sum(genos[i,filled]!=genos[j,filled])
        }else {
          distmat[i,j]=-1
        }
      }
  }

  

  distmat <- as.dist(distmat)
  net <- makeNetSingle(distmat,meta)
  #printGraph(net,"Greens",tab)
  
distmatAllVars <- distmat
netAllVars <- net
  #target:
  #targetDist <- as.dist(sym(read.table("Thies_all_manual.PASS.Cls.miss0.5.LMRG.HAP.vcf.dist.tab.txt",sep="\t")[c1,c1]))
  #printGraph(makeNetSingle(targetDist,meta),"Greens","test")
```








```r
meta <- read.table("Thies_metadata_1701.txt",sep="\t",header=T)
colnames(meta)[1]<-"name"
meta <- meta[meta$name %in% inds,]
row.names(meta) <- inds
meta$Date <- as.Date(as.character(meta$Date),"%d/%m/%Y",origin = "2000-01-01")

net <- netAllVars
net$type = "core"
net$boot = 0
```




```r
#nb: lat=Y, lon=X
meta$GPS <- as.character(meta$GPS)
meanLat <- mean(as.numeric(t(as.data.frame(strsplit(meta$GPS,", ")[3:10]))[,1]))
meanLon <- mean(as.numeric(t(as.data.frame(strsplit(meta$GPS,", ")[3:10]))[,2]))
maxLat <- max(as.numeric(t(as.data.frame(strsplit(meta$GPS,", ")[3:10]))[,1]))
medianLon <- median(as.numeric(t(as.data.frame(strsplit(meta$GPS,", ")[3:10]))[,2]))
#if this isn't Thies you're in trouble...
meta$GPS[1] <- paste(maxLat+0.03,medianLon,sep=", ")
meta$GPS[2] <- paste(maxLat+0.02,medianLon,sep=", ")
gpspos <- matrix(as.numeric(t(as.data.frame(strsplit(meta$GPS,", ")))),ncol = 2)
#reverse for plotting
gpspos <- gpspos[,c(2,1)]
gpscurves <- rep(0,14)
#gpscurves[c(2,4,5)] <- c(-1,-1,20)
gpscurves[c(2,5)] <- c(+0.2,0)
plot.new()
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10-1.png)

```r
# plot(ig,layout=gpspos,vertex.size=5,
# edge.curved=gpscurves,edge.width=as.numeric(factor(net$type,levels=c("boot","core"))[2:15]),
# edge.label.cex=1,edge.label.family="Arial",vertex.label.family="Arial",edge.label=NA)
# #get graph with single-weighted edges

igW <- igA
E(igW)$weight <- 1
gdists <- distances(igW)
distdist <- data.frame(from=character(),
to=character(),
graphdist=integer(),
geodist=numeric(),
stringsAsFactors=F)
for (i in c(2:10)) {
for (j in c(i:10)) {
if (i != j) {
di <- dim(distdist)[[1]]+1
dist = sqrt(abs(gpspos[i,1]-gpspos[j,1])**2+abs(gpspos[i,2]-gpspos[j,2])**2)
write(paste(inds[i],inds[j],gdists[i,j],dist),stderr())
distdist[di,c(1:2)] <- c(inds[i],inds[j])
distdist[di,c(3:4)] <- c(gdists[i,j],dist)
}
}
}
```


```r
opts_chunk$set(fig.width=16, fig.height=16)
opts_chunk$set(dev=c('png','postscript'))
```





```r
ggmap(thies)
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12-1.png)

```r
#igApure <- igA
igA <- as.igraph(netAllVars)
E(igA)$latend <- gpspos[2:10,2]
E(igA)$lonend <- gpspos[2:10,1]
E(igA)$lat <- gpspos[na.omit(netAllVars$ances),2]
E(igA)$lon <- gpspos[na.omit(netAllVars$ances),1]

igf <- fortify(netAllVars)
igf$vertex.names <- as.character(igf$name)
igf$vertex.names[1:10] <- inds
igf$lat <- gpspos[,2]
igf$lon <- gpspos[,1]
igf$latend[!is.na(igf$weight)] <- gpspos[na.omit(netAllVars$ances),2]
igf$lonend[!is.na(igf$weight)] <- gpspos[na.omit(netAllVars$ances),1]
igf
```

```
##          id ances weight       date ances.date year     name vertex.names
## Th086.07  1    NA     NA 0007-10-11       <NA> 2007 Th086.07     Th086.07
## Th106.09  2     1    576 0009-10-26 0007-10-11 2009 Th106.09     Th106.09
## Th106.11  3     2    358 0011-10-27 0009-10-26 2011 Th106.11     Th106.11
## Th117.11  4     3    111 0011-11-02 0011-10-27 2011 Th117.11     Th117.11
## Th132.11  5     2    121 0011-11-11 0009-10-26 2011 Th132.11     Th132.11
## Th134.11  6     3    105 0011-11-14 0011-10-27 2011 Th134.11     Th134.11
## Th162.12  7     2    125 0012-10-10 0009-10-26 2012 Th162.12     Th162.12
## Th196.12  8     5    139 0012-10-25 0011-11-11 2012 Th196.12     Th196.12
## Th230.12  9     5    106 0012-11-09 0011-11-11 2012 Th230.12     Th230.12
## Th074.13 10     2     88 0013-09-23 0009-10-26 2013 Th074.13     Th074.13
##               lat       lon   latend    lonend
## Th086.07 14.83849 -16.92289       NA        NA
## Th106.09 14.82849 -16.92289 14.83849 -16.92289
## Th106.11 14.80636 -16.90927 14.82849 -16.92289
## Th117.11 14.79721 -16.92927 14.80636 -16.90927
## Th132.11 14.80849 -16.91518 14.82849 -16.92289
## Th134.11 14.80596 -16.91990 14.80636 -16.90927
## Th162.12 14.80681 -16.92588 14.82849 -16.92289
## Th196.12 14.79876 -16.93862 14.80849 -16.91518
## Th230.12 14.76843 -16.90779 14.80849 -16.91518
## Th074.13 14.77348 -17.06051 14.82849 -16.92289
```

```r
#ggplot(data = igf,aes(x = lat, y = lon, xend = latend, yend = lonend)) +
library(ggmap)
#thies <- get_map(location = c(lon = meanLon, lat = meanLat), zoom = 12, maptype = 'terrain')
b=0.0

thies <- get_map(location=c(min(gpspos[,1])-b,min(gpspos[,2])-b,max(gpspos[,1])-b,max(gpspos[,2]-b)),zoom=12,maptype = "roadmap")
```

```
## Warning: bounding box given to google - spatial extent only approximate.
```

```
## converting bounding box to center/zoom specification. (experimental)
```

```
## Map from URL : http://maps.googleapis.com/maps/api/staticmap?center=14.80346,-16.984151&zoom=12&size=640x640&scale=2&maptype=roadmap&language=en-EN&sensor=false
```

```r
ggmap(thies)
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12-2.png)

```r
ggmap(thies)+
geom_point(data=igf, aes(x = lon, y = lat), size=5) +
geom_text(data=igf,aes(label = name), size=8) +
geom_edges(data = igf[!is.na(igf$weight),],aes(x = lonend, y = latend, xend = lon, yend = lat),arrow=arrow(length = unit(8, "pt"), type = "closed"),colour="black",size=2) +
theme_blank()
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12-3.png)

```r
igf[,c("lat","lon","latend","lonend")] <- round((igf[,c("lat","lon","latend","lonend")]*2),3)/2
thies <- get_map(location=c(min(gpspos[,1])-b,min(gpspos[,2])-b,max(gpspos[,1])-b,max(gpspos[,2]-b)),zoom=12,maptype = "terrain")
```

```
## Warning: bounding box given to google - spatial extent only approximate.
```

```
## converting bounding box to center/zoom specification. (experimental)
```

```
## Map from URL : http://maps.googleapis.com/maps/api/staticmap?center=14.80346,-16.984151&zoom=12&size=640x640&scale=2&maptype=terrain&language=en-EN&sensor=false
```

```r
ggmap(thies)+
geom_point(data=igf, aes(x = lon, y = lat), size=5) +
geom_text(data=igf,aes(label = name), size=8) +
geom_edges(data = igf[!is.na(igf$weight),],aes(x = lonend, y = latend, xend = lon, yend = lat),arrow=arrow(length = unit(8, "pt"), type = "closed"),colour="black",size=2) +
theme_blank()
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12-4.png)

```r
igf[,c("lat","lon","latend","lonend")] <- round(igf[,c("lat","lon","latend","lonend")],2)

ggmap(thies)+
geom_point(data=igf, aes(x = lon, y = lat), size=5) +
geom_text(data=igf,aes(label = name), size=8) +
geom_edges(data = igf[!is.na(igf$weight),],aes(x = lonend, y = latend, xend = lon, yend = lat),arrow=arrow(length = unit(8, "pt"), type = "closed"),colour="black",size=2) +
theme_blank()
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12-5.png)

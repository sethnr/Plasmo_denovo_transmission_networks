
```r
library(ggplot2)
library(VennDiagram)

opts_chunk$set(fig.width=9, fig.height=9)
opts_chunk$set(dev=c('png'))
```

INDEL counts, disco/haplo100/haplo250

```r
counts <- read.table("vennCounts.txt",col.names=c("type","cat","count"),colClasses = c("character","character","integer"))
icnt <- counts$count[counts$type=="INDEL"]
names(icnt) <- counts$cat[counts$type=="INDEL"]

a1 = icnt['100']+icnt['110']+icnt['111']+icnt['101']
a2 = icnt['010']+icnt['110']+icnt['111']+icnt['011']
a3 = icnt['001']+icnt['011']+icnt['111']+icnt['101']

#draw.triple.venn(area1 = icnt['100'], area2 = icnt['010'], area3 = icnt['001'], 
draw.triple.venn(area1 = a1, area2 = a2, area3 = a3, 
                 n12 = icnt['110']+icnt['111'], n23 = icnt['011']+icnt['111'], n13 = icnt['101']+icnt['111'], 
                 n123 = icnt['111'], category = c("Discovar", "Haplo100bp", "Haplo250bp"), lty = "blank", 
                 fill = c("skyblue", "pink1", "mediumorchid"),euler.d=T,cex=2)
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2-1.png)

```
## (polygon[GRID.polygon.203], polygon[GRID.polygon.204], polygon[GRID.polygon.205], polygon[GRID.polygon.206], polygon[GRID.polygon.207], polygon[GRID.polygon.208], text[GRID.text.209], text[GRID.text.210], text[GRID.text.211], text[GRID.text.212], text[GRID.text.213], text[GRID.text.214], text[GRID.text.215], text[GRID.text.216], text[GRID.text.217], text[GRID.text.218])
```

SNP counts, disco/haplo100/haplo250

```r
scnt <- counts$count[counts$type=="SNP"]
names(scnt) <- counts$cat[counts$type=="SNP"]

a1 = scnt['100']+scnt['110']+scnt['111']+scnt['101']
a2 = scnt['010']+scnt['110']+scnt['111']+scnt['011']
a3 = scnt['001']+scnt['011']+scnt['111']+scnt['101']

#draw.triple.venn(area1 = scnt['100'], area2 = scnt['010'], area3 = scnt['001'], 
draw.triple.venn(area1 = a1, area2 = a2, area3 = a3, 
                 n12 = scnt['110']+scnt['111'], n23 = scnt['011']+scnt['111'], n13 = scnt['101']+scnt['111'], 
                 n123 = scnt['111'], category = c("Discovar", "Haplo100bp", "Haplo250bp"), lty = "blank", 
                 fill = c("skyblue", "pink1", "mediumorchid"),euler.d=T,cex=2)
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-1.png)

```
## (polygon[GRID.polygon.219], polygon[GRID.polygon.220], polygon[GRID.polygon.221], polygon[GRID.polygon.222], polygon[GRID.polygon.223], polygon[GRID.polygon.224], text[GRID.text.225], text[GRID.text.226], text[GRID.text.227], text[GRID.text.228], text[GRID.text.229], text[GRID.text.230], text[GRID.text.231], text[GRID.text.232], text[GRID.text.233], text[GRID.text.234])
```

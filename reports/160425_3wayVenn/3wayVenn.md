
```r
library(ggplot2)
library(VennDiagram)

opts_chunk$set(fig.width=9, fig.height=9)
opts_chunk$set(dev=c('png'))
```



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
                 fill = c("skyblue", "pink1", "mediumorchid"))
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2-1.png)

```
## (polygon[GRID.polygon.43], polygon[GRID.polygon.44], polygon[GRID.polygon.45], polygon[GRID.polygon.46], polygon[GRID.polygon.47], polygon[GRID.polygon.48], text[GRID.text.49], text[GRID.text.50], text[GRID.text.51], text[GRID.text.52], text[GRID.text.53], text[GRID.text.54], text[GRID.text.55], text[GRID.text.56], text[GRID.text.57], text[GRID.text.58])
```


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
                 fill = c("skyblue", "pink1", "mediumorchid"))
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-1.png)

```
## (polygon[GRID.polygon.59], polygon[GRID.polygon.60], polygon[GRID.polygon.61], polygon[GRID.polygon.62], polygon[GRID.polygon.63], polygon[GRID.polygon.64], text[GRID.text.65], text[GRID.text.66], text[GRID.text.67], text[GRID.text.68], text[GRID.text.69], text[GRID.text.70], text[GRID.text.71], text[GRID.text.72], text[GRID.text.73], text[GRID.text.74])
```

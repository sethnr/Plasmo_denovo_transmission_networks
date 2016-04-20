#seqtrack tutorial  

```r
library("phyloTop")
```

```
## Error in library("phyloTop"): there is no package called 'phyloTop'
```

```r
library("phylobase")
```

```
## Error in library("phylobase"): there is no package called 'phylobase'
```

```r
library(knitr)
library(igraph)
library("RColorBrewer")

opts_chunk$set(fig.width=9, fig.height=9)
opts_chunk$set(dev=c('png'))
```



```r
#thies <- read.nexus("sum_ALL.target.2.nexus")
tree <- read.newick("sum_ALL.target.newick")
```

```
## Error in eval(expr, envir, enclos): could not find function "read.newick"
```

```r
tree$tip.label<- gsub("'","",tree$tip.label)
```

```
## Error in tree$tip.label: object of type 'closure' is not subsettable
```

```r
is.rooted(tree)
```

```
## Error in eval(expr, envir, enclos): could not find function "is.rooted"
```

```r
red <- c("Th061.13", "Th095.13", "Th068.12")
blue <- c("Th166.12", "Th245.13", "Th211.13", "Th246.13", "Th092.13")
green <- c("Th230.12","Th196.12","Th106.09","Th074.13","Th106.11","Th117.11","Th134.11","Th086.07","Th162.12","Th132.11")
  
  
tree
```

```
## function (...) 
## constructor_spec(make_tree, ...)
## <environment: namespace:igraph>
```

```r
splitTop(tree,1)
```

```
## Error in eval(expr, envir, enclos): could not find function "splitTop"
```

```r
splitTop(tree,2)
```

```
## Error in eval(expr, envir, enclos): could not find function "splitTop"
```

```r
redtree <- extract.clade(tree,node=34)
```

```
## Error in eval(expr, envir, enclos): could not find function "extract.clade"
```

```r
bluetree <- extract.clade(tree,node=21)
```

```
## Error in eval(expr, envir, enclos): could not find function "extract.clade"
```

```r
greentree <- extract.clade(tree,node=25)
```

```
## Error in eval(expr, envir, enclos): could not find function "extract.clade"
```

```r
# subtree <- greentree
# sackin.phylo(subtree)
# widths(subtree)
# avgLadder(subtree)
# ILnumber(subtree)
# colless.phylo(subtree)
# ladderSizes(subtree)
# ladderShow(subtree)

phyloTop(list(bluetree,greentree,redtree))
```

```
## Error in eval(expr, envir, enclos): could not find function "phyloTop"
```

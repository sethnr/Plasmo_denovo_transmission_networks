library(ggplot2)

rm(CNVs)
rm(CNV_cf)
#value index
V=5
for (file in list.files("./","*CNVs$")) {
  # write(file,stderr())
  tab <- read.table(file,col.names=c("chr","st","en","CN","CNV"))
  set=gsub("_CNVs","",file)
  tab$set=set
  if (!exists("CNVs")) {
    
    CNVs <- tab
    colnames(tab)[V] <- set
    CNVcf <- tab[,c(1:3,V)]    
  } else {
    CNVs <- rbind(CNVs,tab)
    colnames(tab)[V] <- set
    CNVcf <- merge(CNVcf,tab[,c(1:3,V)],by=c("chr","st","en"),all=T)
  }
}
head(CNVcf)
order <- colnames(CNVcf)[4:dim(CNVcf)[[2]]]
order <- sort(order)
CNVcf <- CNVcf[order(as.character(CNVcf$chr)),c("chr","st","en",order)]
write.table(CNVcf,"all_CNVs_merge.txt",sep="\t",row.names=F,col.names=T,quote=F)               



CNVs$set <- factor(CNVs$set,levels=order)
CNVs$chr <- factor(CNVs$chr,levels=sort(levels(CNVs$chr)))

#cnvsCF <- melt(CNVs[,c("chr","st","en","set","CNV")],id.vars=c("chr","st","en","set"))

#ggplot(CNVs,aes(xmin=st,xmax=en,ymin=set,ymax=set,colour=CNV)) + geom_rect(size=1) + geom_rect(data=coronin,ymin=0,ymax=8,colour="black",fill=NA) + facet_grid(chr ~ .) 

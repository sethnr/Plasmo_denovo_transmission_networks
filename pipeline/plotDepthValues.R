library(ggplot2)
library(reshape2)


if(!interactive()) {
  cmd_args = commandArgs(trailingOnly = T);
  prefix = cmd_args[1]
#  depthfile = cmd_args[1]
#  headers = cmd_args[2]
  depthfile = paste(prefix,"depth",sep=".")
  headers = paste(prefix,"samples",sep=".")
} else {  
  depthfile <- "./disco_greenhouse/3D7DD2_greenhouse.depth"
  headers <-  "./disco_greenhouse/3D7DD2_greenhouse.samples"
  prefix <- "3D7DD2_greenhouse"
}

names <- c("CHROM","POS",t(read.table(headers,stringsAsFactors = F)[1,]))

#write(paste(names),stderr())
depths <- read.table(depthfile,sep="\t",col.names = names)
#write.table(head(depths),stderr())
depths.m <- melt(depths,id.vars = c("CHROM","POS"))

depthplot <- ggplot(depths.m,aes(x=POS,y=value,group=variable,colour=variable)) + 
  geom_line() + facet_grid(variable ~ CHROM,scales="free_x")

W=20
H=1
ggsave(paste(prefix,".png",sep=""),depthplot,width =W,height=H*length(levels(depths.m$variable)))

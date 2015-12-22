library(ggplot2)
library(reshape2)
library(zoo)

popSizes <- t(read.table("Asexuals_100.txt",sep="\t"))
popSizes <- as.data.frame(popSizes)
popSizes$t <- c(1:dim(popSizes)[[1]])

popSizes.m <- melt(popSizes,id.vars=c("t"))
ggplot(popSizes.m,aes(x=t,y=value,group=variable)) + geom_line(alpha=0.1) + scale_y_log10()


popSizes$mean <- apply(popSizes[,1:100],1,FUN=mean,na.rm=T)
popSizes$max <- apply(popSizes[,1:100],1,FUN=max,na.rm=T)

ggplot(popSizes,aes(x=t,y=mean)) + geom_line() + geom_line(aes(x=t,y=max),colour="red") + scale_y_log10()

smoothmax <- loess(value ~ t, span=0.05, popSizes.m)
popSizes$smooth <- predict(smoothmax,data.frame(t=popSizes$t))

#loess smooth max values, repredict from dates
#popSizes$smooth*5e6
m=5e6
ggplot(popSizes,aes(x=t,y=smooth*m)) + geom_line() + geom_line(aes(x=t,y=max*m),colour="red") + scale_y_log10()
maxtab <- cbind(popSizes$t,predict(smoothmax,data.frame(t=popSizes$t))*m)
colnames(maxtab) <- c("T","C")

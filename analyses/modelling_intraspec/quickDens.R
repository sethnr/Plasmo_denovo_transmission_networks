afs_small <- data.frame(t(read.table("Scratch/gits/pfdisco/analyses/modelling_intraspec/qstages.C1.5m.L0.01mb.G36.mtp.AFs.txt",sep="\t",row.names=1)))
colnames(afs_small)<-c("mth","halfyr","yr")
afs_large <- data.frame(t(read.table("Scratch/gits/pfdisco/analyses/modelling_intraspec/qstages.C100.0m.L0.01mb.G36.mtp.AFs.txt",sep="\t",row.names=1)))
colnames(afs_large)<-c("mth","halfyr","yr")
plot(density(afs_large$mth),col="red") + lines(density(afs_small$yr),col="blue")


plot(density(subset(afs_large,mth>0.01)$mth),col="red",xlim=c(0,0.1)) + 
  lines(density(subset(afs_small,halfyr>0.01)$halfyr),col="blue") +
  lines(density(subset(afs_small,yr>0.01)$yr),col="dark blue")


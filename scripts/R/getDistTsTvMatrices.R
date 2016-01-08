library(ggplot2)
library(ape)
library(stringdist)
library(reshape2)


fileroot <- commandArgs(TRUE)[1]

ped=paste(fileroot,'ped',sep='.')
outpng=paste(fileroot,'png',sep='.')
outnex=paste(fileroot,'nexus',sep='.')
outtab_dist=paste(fileroot,'dist.tab.txt',sep='.')
outtab_tstv=paste(fileroot,'tstv.tab.txt',sep='.')
outtab_id=paste(fileroot,'in-del.tab.txt',sep='.')
outtab_IS=paste(fileroot,'indel-snp.tab.txt',sep='.')
outtab_taat=paste(fileroot,'taat.tab.txt',sep='.')

#read in ped file and get hamming distances:
genos <- read.table(ped,colClasses="character")
inds <- genos[,1]
genos <- genos[,seq(7,dim(genos)[[2]],2)]
rownames(genos)=inds

distmat = matrix(nrow=length(inds),ncol=length(inds))
colnames(distmat) = inds
rownames(distmat) = inds

#write("calculating distance matrix",stderr())
#for (i in rownames(genos)){
#    for (j in rownames(genos)){
#    filled = intersect(which(genos[i,] != 0), which(genos[j,] !=0))
#    write(paste("calculating",i,"v",j),stderr())
#    write(length(filled),stderr())
#    distmat[i,j] = sum(genos[i,filled]!=genos[j,filled])
#    }
#}

#write.table(distmat,stderr())
#distmat[lower.tri(distmat)] <- ''
#write.table(distmat,sep="\t",file=outtab_dist,quote=F)


tstvmat = matrix(nrow=length(inds),ncol=length(inds))
colnames(tstvmat) = inds
rownames(tstvmat) = inds

idmat = matrix(nrow=length(inds),ncol=length(inds))
colnames(idmat) = inds
rownames(idmat) = inds

ISmat = matrix(nrow=length(inds),ncol=length(inds))
colnames(ISmat) = inds
rownames(ISmat) = inds

TAATmat = matrix(nrow=length(inds),ncol=length(inds))
colnames(TAATmat) = inds
rownames(TAATmat) = inds


snpClasses <- expand.grid(c("A","C","T","G"),c("A","C","T","G"))
snpCombs <- paste(snpClasses$Var1,snpClasses$Var2,sep="/")
snpCombs <- snpCombs[snpClasses$Var1 != snpClasses$Var2]
snpClasses <- rep("tv",length(snpCombs))
names(snpClasses) <- snpCombs
snpClasses["A/G"] <- "ts"
snpClasses["G/A"] <- "ts"
snpClasses["C/T"] <- "ts"
snpClasses["T/C"] <- "ts"


write("calculating TvTs matrix",stderr())
for (i in rownames(genos)){
    for (j in rownames(genos)){
        filled = intersect(which(genos[i,] != 0), which(genos[j,] !=0))
#    	write(paste("calculating",i,"v",j),stderr())
#      	write(length(filled),stderr())
	
	distmat[i,j] = sum(genos[i,filled]!=genos[j,filled])
	
	#CALC Ts:Tv ratio and IN:DEL ratio
	ts<-0
    	tv<-0
	taat<-0
    	IN <- 0
    	DEL <- 0
    	na=0
    	for (f in filled) {
    	    if ( genos[i,f] != genos[j,f] ) {
#	        write(paste(	genos[i,f],genos[j,f], 
#	  		genos[i,f] != genos[j,f],
#			#paste(genos[i,f],genos[j,f],sep="/"),
#			snpClasses[paste(genos[i,f],genos[j,f],sep="/")],
#			sep="\t"),stderr())
		varID <- paste(genos[i,f],genos[j,f],sep="/")

	  	tstv <- snpClasses[varID]
	  	if (is.na(tstv)) { 
 	     	   if(nchar(genos[i,f]) > nchar(genos[j,f])) {DEL <- DEL+1
	     	   } else if(nchar(genos[i,f]) < nchar(genos[j,f])) {IN <- IN+1
	     	   } else {na <- na+1}
		} else if ( snpClasses[varID] == "tv" ) {
		    tv <- tv+1
	  	} else if ( snpClasses[varID] == "ts" ) {
	    	    ts <- ts+1}
		if (varID %in% c("T/A","A/T")) {taat <- taat+1}
	    }
	}
    
	write(paste(i,j,distmat[i,j],ts,tv,IN,DEL,na,sep="\t"),stderr())
    	tstvmat[i,j] = round(tv/ts,2)
    	idmat[i,j] = round(IN/DEL,2)
	ISmat[i,j] = round((tv+ts)/(IN+DEL),2)
    	TAATmat[i,j] = round(taat/(tv+ts),2)
    }
}

distmat[lower.tri(distmat)] <- ''
write.table(distmat,sep="\t",file=outtab_dist,quote=F)
write.table(tstvmat,sep="\t",file=outtab_tstv,quote=F)
write.table(idmat,sep="\t",file=outtab_id,quote=F)
write.table(ISmat,sep="\t",file=outtab_IS,quote=F)
write.table(TAATmat,sep="\t",file=outtab_taat,quote=F)




# idmat = matrix(nrow=length(inds),ncol=length(inds))
# colnames(idmat) = inds
# rownames(idmat) = inds

# write("calculating IN:DEL matrix",stderr())
# for (i in rownames(genos)){
#     for (j in rownames(genos)){
#     filled = intersect(which(genos[i,] != 0), which(genos[j,] !=0))
#     write(paste("calculating",i,"v",j),stderr())
#     write(length(filled),stderr())
#     IN<-0
#     DEL<-0
#     for (f in filled) {
#     	if (genos[i,f]==genos[j,f]){ next
#	}else {
#	  if(nchar(genos[i,f]) > nchar(genos[j,f])) {DEL +=1}
#	  if(nchar(genos[i,f]) < nchar(genos[j,f])) {IN +=1}
#	}
#    }
#    idmat[i,j] = IN/DEL
#    }
#}

#write.table(idmat,sep="\t",file=outtab_id,quote=F)


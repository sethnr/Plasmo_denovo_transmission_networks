
SAMPLE=$1
OLDNAME=$2

cd $SAMPLE

for GZVCF in *vcf.gz
do
  gunzip $GZVCF
  perl -i -pe "s/${OLDNAME}$/${SAMPLE}/gi" ${GZVCF/.gz/}
  head -n 40 ${GZVCF/.gz/} | grep \#CHROM
  bgzip ${GZVCF/.gz/}
  tabix -pvcf ${GZVCF}
done
  

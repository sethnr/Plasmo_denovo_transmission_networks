#!/bin/bash

NUCMER=$1
DISCO=$2
NUCMERBED=$3
#DISCOBED=$4
REGION=$4
OUTFILE=$5

echo "setup"
echo "F1" $NUCMER
echo "F2" $DISCO
echo "bed" $NUCMERBED
echo "region" $REGION
echo "outfile" $OUTFILE

#get abspaths
NUCMER=`readlink -f  $NUCMER`
DISCO=`readlink -f  $DISCO`
NUCMERBED=`readlink -f  $NUCMERBED`
#DISCOBED=`readlink -f  $DISCOBED`

echo "F1" $NUCMER
echo "F2" $DISCO
echo "bed" $NUCMERBED
echo "region" $REGION
echo "outfile" $OUTFILE


mkdir $OUTFILE
cd $OUTFILE

#make local links
ln -s $NUCMER `basename $NUCMER`
NUCMER=`basename $NUCMER`
ln -s $DISCO `basename $DISCO`
DISCO=`basename $DISCO`
ln -s $NUCMERBED `basename $NUCMERBED`
NUCMERBED=`basename $NUCMERBED`
#ln -s $DISCOBED `basename $DISCOBED`
#DISCOBED=`basename $DISCOBED`


echo "parsing matching regions"
# vcftools --gzvcf $NUCMER --recode --out ${NUCMER/.vcf*/}.cf --bed $DISCOBED
# echo "done"
# mv ${NUCMER/.vcf*/}.cf.recode.vcf ${NUCMER/.vcf*/}.cf.vcf
# vcftools --gzvcf $DISCO --recode --out ${DISCO/.vcf*/}.cf --bed $NUCMERBED
# mv ${DISCO/.vcf*/}.cf.recode.vcf ${DISCO/.vcf*/}.cf.vcf

tabix -pvcf $DISCO
tabix -pvcf $NUCMER

tabix -h $DISCO $REGION > ${DISCO/.vcf*/}.cf.vcf
tabix -h $NUCMER $REGION > ${NUCMER/.vcf*/}.cf.vcf

echo "comparing vcfs"
echo python $DISCO1/scripts/cfs/cfVCFs.py \
  --v1 ${NUCMER/.vcf*/}.cf.vcf \
  --v2 ${DISCO/.vcf*/}.cf.vcf \
  -o ${OUTFILE}.vcf
python $DISCO1/scripts/cfs/cfVCFs.py \
  --v1 ${NUCMER/.vcf*/}.cf.vcf \
  --v2 ${DISCO/.vcf*/}.cf.vcf \
  -o ${OUTFILE}.vcf  1>/dev/null


echo "adding STR information"
echo python $DISCO1/scripts/cfs/addSTRsToVCF.py -v ${OUTFILE}.vcf -s $DISCO1/analyses/STR_calling/STRs_mreps_Pf3D7_v3.R0.txt
python $DISCO1/scripts/cfs/addSTRsToVCF.py -v ${OUTFILE}.vcf -s $DISCO1/analyses/STR_calling/STRs_mreps_Pf3D7_v3.R0.txt
# mv ${OUTFILE}.STRs.vcf ${OUTFILE}.vcf



echo "calculating match statistics"
python $DISCO1/scripts/cfs/getVcfMatchStats.py -v ${OUTFILE}.STRs.vcf > ${OUTFILE}.txt


cp ${OUTFILE}.txt ../${OUTFILE}.txt
cp ${OUTFILE}.STRs.vcf ../${OUTFILE}.vcf

exit 0

#python ../scripts/getDepthHet.py -v fakev3D7_all.inDD2.recode.cfout.vcf.gz -d fakev3D7.depth.gz -r Pf3D7_07_v3:450000-750000 -b10 > depth_v_het_10kb_150629.txt



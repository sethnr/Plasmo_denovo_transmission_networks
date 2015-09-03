#!/bin/bash

OUT=$1
shift

#SAMPLE_NO=`ls $OUT | wc -l`
#VCF_NO=`find $OUT | grep vcf.gz$ | wc -l`

#TOMERGE=`find $OUT | grep vcf.gz$ `


#check if it is already there (i.e. SNPs already calculated)
ls ${OUT}_all.vcf.gz
rc=$?;
if [[ $rc == 0 ]];
then
    echo "MERGE/CONCAT OUTPUT ALREADY PRESENT, NOT RUNNING";
    exit 0;
fi


#if [[ $SAMPLE_NO != $VCF_NO ]]
#   then
cd $OUT
SAMPLES=`ls -d */`
cd ../

for SAMPLE in $SAMPLES
do
    SAMPLE=${SAMPLE/\/}
    if [[ -d "${OUT}/${SAMPLE}" ]]
    then
	echo "concatting"
	ls ${OUT}/${SAMPLE}/*vcf.gz | head -n 10
	echo '...'
	echo "to ${OUT}/${SAMPLE}_concat.vcf.gz"
	
	vcf-concat -s 1 ${OUT}/${SAMPLE}/*vcf.gz | vcf-sort > ${OUT}/${SAMPLE}_concat.vcf
	bgzip -f ${OUT}/${SAMPLE}_concat.vcf
	tabix -pvcf ${OUT}/${SAMPLE}_concat.vcf.gz
    fi
done
TOMERGE=`ls ${OUT}/*concat.vcf.gz`
#fi



echo "merging"
echo $TOMERGE
echo "to ${OUT}_all.vcf"
vcf-merge -d -s $TOMERGE > ${OUT}_all.vcf
bgzip ${OUT}_all.vcf
tabix -pvcf ${OUT}_all.vcf.gz

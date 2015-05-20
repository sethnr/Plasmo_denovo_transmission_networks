#!/bin/bash

OUT=$1
shift

SAMPLE_NO=`ls $OUT | wc -l`
VCF_NO=`find $OUT | grep vcf.gz$ | wc -l`

TOMERGE=`find $OUT | grep vcf.gz$ `

if [[SAMPLE_NO == VCF_NO]]
   then
       for SAMPLE in `ls $OUT `
       do
	   echo "concatting"
	   ls ${OUT}/${SAMPLE}/*/*vcf.gz
	   echo "to" ${OUT}"/"${SAMPLE}"_concat.vcf.gz"

	   vcf-concat ${OUT}/${SAMPLE}/*/*vcf.gz > ${OUT}/${SAMPLE}_concat.vcf
	   bgzip ${OUT}/${SAMPLE}_concat.vcf
	   tabix -pvcf ${OUT}/${SAMPLE}_concat.vcf.gz
       done
   TOMERGE=`ls ${OUT}/*concat.vcf.gz`
fi

echo "merging"
echo $TOMERGE
echo "to ${OUT}_all.vcf"
vcf-merge -s $TOMERGE > ${OUT}_all.vcf
bgzip ${OUT}_all.vcf
tabix -pvcf ${OUT}_all.vcf.gz

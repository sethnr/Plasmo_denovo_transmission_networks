#!/bin/bash

VCF=$1
SAMPLES=$2
FASTA=$3


echo python $DISCO1/scripts/vcf/mergeLanesVCF.py -v $VCF
#python $DISCO1/scripts/vcf/mergeLanesVCF.py -v $VCF

VCF=${VCF/.gz/}

echo python $DISCO1/scripts/vcf/replaceNamesVCF.py -v ${VCF/.vcf/.LMRG.vcf} -f $SAMPLES
python $DISCO1/scripts/vcf/replaceNamesVCF.py -v ${VCF/.vcf/.LMRG.vcf} -f $SAMPLES

echo $DISCO1/scripts/splitVcfBySample.sh ${VCF/.vcf/.LMRG.RENAME.vcf} $FASTA
$DISCO1/scripts/splitVcfBySample.sh ${VCF/.vcf/.LMRG.RENAME.vcf} $FASTA

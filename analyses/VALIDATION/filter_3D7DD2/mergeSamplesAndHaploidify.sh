#!/bin/bash

VCF=$1
SAMPLES=$2
OUTFILE=$3
#FASTA=$3

SAMPLES_ABS=`readlink -f ${SAMPLES}`
SAMPLES=`basename $SAMPLES`

VCF_ABS=`readlink -f ${VCF}`
VCF=`basename $VCF`



mkdir tmp_${OUTFILE}_mergeHap
cd tmp_${OUTFILE}_mergeHap

ln -s $SAMPLES_ABS
ln -s $VCF_ABS

echo python $DISCO1/scripts/vcf/mergeLanesVCF.py -v $VCF
python $DISCO1/scripts/vcf/mergeLanesVCF.py -v $VCF

VCF=${VCF/.gz/}

echo python $DISCO1/scripts/vcf/replaceNamesVCF.py -v ${VCF/.vcf/.LMRG.vcf} -f $SAMPLES
python $DISCO1/scripts/vcf/replaceNamesVCF.py -v ${VCF/.vcf/.LMRG.vcf} -f $SAMPLES

echo python $DISCO1/scripts/vcf/haploidifyVCF.py -v ${VCF/.vcf/.LMRG.RENAME.vcf}
python $DISCO1/scripts/vcf/haploidifyVCF.py -v ${VCF/.vcf/.LMRG.RENAME.vcf}

cp ${VCF/.vcf/.LMRG.RENAME.HAPLOID.vcf} ../${VCF/.vcf/.LMRG.HAP.vcf}

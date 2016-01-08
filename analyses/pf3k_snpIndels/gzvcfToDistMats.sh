#!/bin/bash

VCF=$1
OUTFILE=$2

VCF_ABS=`readlink -f ${VCF}`
VCF=`basename $VCF`

if [[ -z $OUTFILE ]] 
then
    OUTFILE=$VCF
fi

PLINK=/broad/software/free/Linux/redhat_6_x86_64/pkgs/plink_1.90b/plink
PEDTOTREE=$DISCO1/scripts/R/getDistTsTvMatrices.R

mkdir tmp_${OUTFILE}_NJnex
cd tmp_${OUTFILE}_NJnex

ln -s $VCF_ABS

vcftools --plink --out ${VCF/.vcf.gz/} --gzvcf $VCF
if [[ $? != 0 ]]; then exit 1 ; fi

VCF=${VCF/.gz/}

Rscript $PEDTOTREE ${VCF/.vcf/}
if [[ $? != 0 ]]; then exit 1 ; fi


cp *.tab.txt ../

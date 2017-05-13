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
PEDTOTREE=$DISCO1/scripts/R/getNJTreePed.R

mkdir tmp_${OUTFILE}_NJnex
cd tmp_${OUTFILE}_NJnex

ln -s $VCF_ABS

vcftools --plink --out ${VCF/.vcf/} --gzvcf $VCF --maf 0.01 --max-maf 0.99
if [[ $? != 0 ]]; then exit 1 ; fi

Rscript $PEDTOTREE ${VCF/.vcf/}
if [[ $? != 0 ]]; then exit 1 ; fi

cp ${VCF/.vcf/.nexus} ../${OUTFILE}.NJ.nexus
cp ${VCF/.vcf/.png} ../${OUTFILE}.NJ.png
cp ${VCF/.vcf/.dist.tab.txt} ../${OUTFILE}.dist.tab.txt

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
GETALTAB=$WORK/idi_broad_scripts/vcf/getAlleleTableVCF.py
PEDTOTREE=$PWD/getParsTreePed.R

mkdir tmp_${OUTFILE}_PARSnex
cd tmp_${OUTFILE}_PARSnex

ln -s $VCF_ABS

#vcftools --plink --out ${VCF/.vcf.gz/} --gzvcf $VCF
python $GETALTAB -M -1 -v $VCF > ${VCF/.vcf.gz/.alleles.tab}

if [[ $? != 0 ]]; then exit 1 ; fi

Rscript --vanilla $PEDTOTREE ${VCF/.vcf.gz/.alleles}
if [[ $? != 0 ]]; then exit 1 ; fi

cp ${VCF/.vcf.gz/.nexus} ../${OUTFILE}.PARS.nexus
cp ${VCF/.vcf.gz/.png} ../${OUTFILE}.PARS.png
cp ${VCF/.vcf.gz/.dist.tab.txt} ../${OUTFILE}.dist.tab.txt

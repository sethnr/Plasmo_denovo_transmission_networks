#!/bin/bash

VCF=$1

vcftools --plink --out ${VCF/.vcf/} --vcf $VCF
Rscript $DISCO1/scripts/R/getDistTsTvMatrices.R ${VCF/.vcf/}

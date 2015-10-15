#!/bin/bash

VCF=$1

$DISCO1/scripts/diploidifyVCF.py -v $VCF
VCF=${VCF/.gz/}

vcftools --012 --vcf ${VCF/.vcf/.DIPLOID.vcf} --out ${VCF/.vcf/}

Rscript $DISCO1/scripts/R/plotHaplotypes.R ${VCF/.vcf/}

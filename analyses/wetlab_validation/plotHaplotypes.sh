#!/bin/bash

VCF=$1

$DISCO1/scripts/haploidifyVCF.py -v $VCF
VCF=${VCF/.gz/}

vcftools --012 --vcf ${VCF/.vcf/.HAPLOID.vcf} --out ${VCF/.vcf/}

#Rscript $DISCO1/scripts/R/plotHaplotypes.R ${VCF/.vcf/.012/} ${VCF/.vcf/.012.indv/}  ${VCF/.vcf/.012.pos/}

#!/bin/bash

VCF=$1

#python $DISCO1/scripts/vcf/haploidifyVCF.py -v $VCF
#VCF=${VCF/.gz/}

#vcftools --012 --vcf ${VCF/.vcf/.DIPLOID.vcf} --out ${VCF/.vcf/}
#python getAlleleTableVCF.py -v ${VCF/.vcf/.HAPLOID.vcf} > ${VCF/.vcf/.alleles.tab.txt}
python getAlleleTableVCF.py -v ${VCF} | sort -u > ${VCF/.vcf/.alleles.tab.txt}

Rscript $DISCO1/scripts/R/plotHaplotypes.R ${VCF/.vcf/.alleles.tab.txt}

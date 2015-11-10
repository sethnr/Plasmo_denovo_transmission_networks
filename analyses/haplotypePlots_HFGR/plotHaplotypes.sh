#!/bin/bash

VCF=$1

#python $DISCO1/scripts/vcf/haploidifyVCF.py -v $VCF
#VCF=${VCF/.gz/}

#python $DISCO1/scripts/vcf/getAlleleTableVCF.py -v ${VCF} | sort -u > ${VCF/.vcf/.alleles.tab.txt}
python $DISCO1/scripts/vcf/getAlleleTableVCF.py -v ${VCF} > ${VCF/.vcf/.alleles.tab.txt}

Rscript $DISCO1/scripts/R/plotHaplotypes.R ${VCF/.vcf/.alleles.tab.txt}

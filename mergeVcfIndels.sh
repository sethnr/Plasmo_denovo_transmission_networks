#!/bin/bash

OUT=$
shift

vcf-merge $@ > merge.vcf
vcftools --keep-only-indels --vcf merge.vcf --recode --out merge_indels
mv merge_indels.recode.vcf merge_indels.vcf
bgzip merge_indels.vcf
tabix -pvcf merge_indels.vcf.gz
tabix -h merge_indels.vcf.gz $FINELOCS >  ${OUT}.vcf
rm merge.vcf \
   merge.vcf.vcfidx \
   merge_indels.log \
   merge_indels.vcf.gz \
   merge_indels.vcf.gz.tbi

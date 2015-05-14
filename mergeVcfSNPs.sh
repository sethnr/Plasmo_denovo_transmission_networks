#!/bin/bash

OUT=$1
shift
LOCS=$1
shift

echo "MERGING VCFS AND PARSING SNPS"
echo "OUT: " $OUT
echo "LOCS: " $LOCS

#vcf-merge $@ 1> merge.vcf 2>/dev/null
./mergeConcatVcfs.sh $OUT

vcftools --remove-indels --vcf ${OUT}_all.vcf --recode --out merge_snps
mv merge_snps.recode.vcf merge_snps.vcf
bgzip merge_snps.vcf
tabix -pvcf merge_snps.vcf.gz


# replace commas with spaces for tabix
LOCSP="${LOCS//\,/\ }"
tabix -h merge_snps.vcf.gz $LOCSP >  ${OUT}_SNPs.vcf


rm merge.vcf \
   merge.vcf.vcfidx \
   merge_snps.log \
   merge_snps.vcf.gz \
   merge_snps.vcf.gz.tbi

exit 0
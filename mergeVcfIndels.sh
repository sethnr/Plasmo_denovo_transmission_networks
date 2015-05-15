#!/bin/bash

OUT=$1
shift
LOCS=$1
shift

echo "MERGING VCFS AND CALLING INDELS"
echo "OUT: " $OUT
echo "LOCS: " $LOCS

SCRIPTDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
# vcf-merge $@ 1> merge.vcf 2>/dev/null
echo ${SCRIPTDIR}/mergeConcatVcfs.sh $OUT

${SCRIPTDIR}/mergeConcatVcfs.sh $OUT

vcftools --keep-only-indels --vcf ${OUT}_all.vcf --recode --out merge_indels
mv merge_indels.recode.vcf merge_indels.vcf
bgzip merge_indels.vcf
tabix -pvcf merge_indels.vcf.gz


# replace commas with spaces for tabix
LOCSP="${LOCS//\,/ }"
tabix -h merge_indels.vcf.gz $LOCSP > ${OUT}_STRs.vcf


rm merge.vcf \
   merge.vcf.vcfidx \
   merge_indels.log \
   merge_indels.vcf.gz \
   merge_indels.vcf.gz.tbi

#UGLY HACK: (discovar leaves region=<nothing> tags in header, STR parser chokes)
perl -i -ne 'print $_ unless $_ =~ m/DiscovarRegion/gi' ${OUT}_STRs.vcf


exit 0

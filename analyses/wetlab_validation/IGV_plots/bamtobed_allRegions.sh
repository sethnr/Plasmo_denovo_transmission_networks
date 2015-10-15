#!/bin/bash

BEDS=""
for BAM in $@;
do 
    echo bedtools bamtobed -i $BAM \> ${BAM/.bam/.bed}
    bedtools bamtobed -i $BAM > ${BAM/.bam/.bed}
    BEDS="${BEDS} ${BAM/.bam/.bed}"
done
cat $BEDS | sort -k1,1 -k2,2g > allRegions.bed
bedtools merge -d 5000 -i allRegions.bed > allRegions.2.bed
mv allRegions.2.bed allRegions.bed

perl -i -lane 'print join("\t",@F[0..2],$F[0]."_".$F[1]."_".$F[2]);' allRegions.bed

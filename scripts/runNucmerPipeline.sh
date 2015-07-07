#!/bin/bash

FASTA1=$1
FASTA2=$2
OUT=$3

echo "running nucmer"
nucmer -c 1000 -o -p $OUT --nooptimize -b 20000 $FASTA1 $FASTA2

echo "parsing SNPs/INDELs"
show-snps -TC ${OUT}.delta > ${OUT}.snps

echo "mummer 2 bed"
#make bed of matching regions
perl $DISCO1/hacks/mummer2bed.pl ${OUT}.coords

echo "mummer 2 VCF"
#make two vcfs from mummer SNP output
python $DISCO1/scripts/mummer2vcf_plus1error.py -m ${OUT}.snps

echo "done"

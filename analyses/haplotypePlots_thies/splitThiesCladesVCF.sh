#!/bin/bash

VCF=$1


VCF_ABS=`readlink -f ${VCF}`
VCF=`basename $VCF`


if [[ -z $OUTFILE ]]
then
    OUTFILE=$VCF
fi

mkdir tmp_${OUTFILE}_clades
cd tmp_${OUTFILE}_clades

ln -s $VCF_ABS

REF="${WORK}/refs/Pf3D7_v3.fasta"
GATK=/humgen/gsa-hpprojects/GATK/bin/GenomeAnalysisTK-3.4-0-g7e26428

java -jar $GATK/GenomeAnalysisTK.jar \
   -T SelectVariants \
   -R $REF \
   -V ${VCF} \
   -o ${VCF/.vcf/.clade3.vcf} \
   -env \
   -sn Th086.07 -sn Th106.09 -sn Th106.11 -sn Th117.11 \
   -sn Th132.11 -sn Th134.11 -sn Th162.12 -sn Th196.12 \
   -sn Th230.12 -sn Th074.13 \
   -trimAlternates

vcftools  --vcf ${VCF/.vcf/.clade3.vcf} --maf 0.1 --max-maf 0.9 \
  --out  ${VCF/.vcf/.clade3} --recode 

java -jar $GATK/GenomeAnalysisTK.jar \
   -T SelectVariants \
   -R $REF \
   -V ${VCF} \
   -o ${VCF/.vcf/.clade2.vcf} \
   -env \
   -sn Th166.12 -sn Th092.13 -sn Th211.13 -sn Th245.13 \
   -sn Th246.13 \
   -trimAlternates

vcftools  --vcf ${VCF/.vcf/.clade2.vcf} --maf 0.1 --max-maf 0.9 \
  --out  ${VCF/.vcf/.clade2} --recode 


java -jar $GATK/GenomeAnalysisTK.jar \
   -T SelectVariants \
   -R $REF \
   -V ${VCF} \
   -o ${VCF/.vcf/.clade1.vcf} \
   -env \
   -sn Th068.12 -sn Th061.13 -sn Th095.13 \
   -trimAlternates

vcftools  --vcf ${VCF/.vcf/.clade1.vcf} --maf 0.1 --max-maf 0.9 \
  --out  ${VCF/.vcf/.clade1} --recode 


#!/bin/bash

VCF=$1
shift

FASTA=$1
shift

GATK=/humgen/gsa-hpprojects/GATK/bin/GenomeAnalysisTK-3.5-0-g36282e4
#if any arguments left after vcf file
if [[ "$#" > 0 ]]
then
    echo "passed samples"
    SAMPLES=$@
else
    echo "parsed samples"
    if [[ $VCF == *".gz" ]]
     then 
 	SAMPLES=`zcat $VCF | head -n 1000 | grep \#CHR | cut -f 10- `
     else
	SAMPLES=`head -n 1000 ${VCF} | grep \#CHR | cut -f 10- `
    fi

    # TO DO: 
    # parse out ref from VCF :
    # file:///seq/plasmodium/sredmond/refs/Pf3D7_v3.fasta
fi 

for SAMPLE in $SAMPLES 
do
echo "SAMPLE = " $SAMPLE
OUT=${VCF/.vcf/.${SAMPLE}.vcf}
echo $VCF "->" $OUT
echo "REF " $FASTA
echo java -jar ${GATK}/GenomeAnalysisTK.jar \
   -R $FASTA \
   -T SelectVariants \
   -V $VCF \
   -o $OUT  \
   -sn $SAMPLE \
   -env \
#   -ploidy 1 \
#   -trimAlternates
java -jar ${GATK}/GenomeAnalysisTK.jar \
   -R $FASTA \
   -T SelectVariants \
   -V $VCF \
   -o $OUT  \
   -sn $SAMPLE \
   -env \
#   -ploidy 1 \
#   -trimAlternates

#   -ef \
done

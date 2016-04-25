#!/bin/bash

while getopts "d:R:l:o:V:" opt; do
  case $opt in
    l) LOCUS=$OPTARG;;
    R) REF=${OPTARG};;
    o) OUTFILE=$OPTARG;;
    V) VCF=$OPTARG;;
   \?) echo "Invalid option: -$OPTARG" >&2 ;;
  esac
done
shift $(expr $OPTIND - 1 )

echo $REF

#echo samtools faidx $REF $LOCUS \> ${LOCUS}.fasta
#samtools faidx $REF $LOCUS > ${LOCUS}.fasta
#tabix -h $VCF $LOCUS > ${LOCUS}.vcf
#bgzip ${LOCUS}.vcf
#tabix -pvcf ${LOCUS}.vcf.gz
#splitVcfBySample.sh ${LOCUS}.vcf.gz $REF
#samtools faidx $REF $LOCUS > locus.fasta
#tabix -h $VCF $LOCUS > locus.vcf
#bgzip locus.vcf
#tabix -pvcf locus.vcf.gz
#splitVcfBySample.sh locus.vcf.gz $REF

if [[ "$#" > 0 ]]
then
    echo "passed samples"
    SAMPLES=$@
else

    if [[ $VCF == *".gz" ]]
    then 
	SAMPLES=`zcat $VCF | head -n 1000 | grep \#CHR | cut -f 10- `
    else
	SAMPLES=`head -n 1000 ${VCF} | grep \#CHR | cut -f 10- `
    fi
fi

rm locus.fasta
for SAMP in $SAMPLES ; 
do  

  echo samtools faidx $REF $LOCUS \| vcf-consensus -s $SAMP $VCF \> ${SAMP}.fasta
  samtools faidx $REF $LOCUS | vcf-consensus -s $SAMP $VCF 1> ${SAMP}.fasta 2>/dev/null
  perl -i -pe "s/\>/\>${SAMP}\t/" ${SAMP}.fasta
  cat ${SAMP}.fasta >> locus.fasta
done

clustalw2 -infile=locus.fasta -outfile=locus.nexus -output=nexus


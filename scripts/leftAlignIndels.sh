#!/bin/bash


PICARD=/seq/software/picard/1.99/bin
GATK=/humgen/gsa-hpprojects/GATK/bin/GenomeAnalysisTK-3.4-0-g7e26428

while getopts "v:f:o:" opt; do
  case $opt in
    v) VCF=${OPTARG};;
    f) REF=$OPTARG;;
    o) OUTFILE=$OPTARG;;
   \?) echo "Invalid option: -$OPTARG" >&2 ;;
  esac
done

if [[ -z "$OUTFILE" ]]
then 
    echo "OUTFILE not found"
     OUTFILE=$VCF
fi

if [[ $VCF == $OUTFILE ]];
then
  echo "modifying in-place"
  java -jar ${GATK}/GenomeAnalysisTK.jar \
   -T LeftAlignAndTrimVariants \
   -R $REF \
   --variant $VCF \
   -o tmp_${OUTFILE}
   mv tmp_${OUTFILE} $OUTFILE
else
  java -jar ${GATK}/GenomeAnalysisTK.jar \
   -T LeftAlignAndTrimVariants \
   -R $REF \
   --variant $VCF \
   -o $OUTFILE
fi

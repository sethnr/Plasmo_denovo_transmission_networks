#!/bin/bash

#NB GATK3.4 requires Java 1.8

LOCIFILE="NULL"
BAMSFILE="NULL"
PLOIDY=1
OUTFILE="GATK_out.vcf"
VARONLY="-out_mode EMIT_ALL_SITES"

PICARD=/seq/software/picard/1.99/bin
GATK=/humgen/gsa-hpprojects/GATK/bin/GenomeAnalysisTK-3.4-0-g7e26428

while getopts "B:b:R:l:p:o:V" opt; do
  case $opt in
    b) BAM=${OPTARG} ;;
    B) BAMSFILE=$OPTARG ;;
#    l) LOCUS=$OPTARG ;;
    l) LOCIFILE=$OPTARG ;;
    p) PLOIDY=$OPTARG;;
    R) REF=$OPTARG;;
    o) OUTFILE=$OPTARG;;
    V) VARONLY="";;
   \?) echo "Invalid option: -$OPTARG" >&2 ;;
  esac
done

#BAMSFILE must end in '.list'
if [[ $BAMSFILE!="NULL" ]]; 
  then 
  BAM=$BAMSFILE
fi

#check ref.dict file is present, otherwise create
if [ ! -f ${REF//fasta/dict} ]
  then
  java -jar ${PICARD}/CreateSequenceDictionary.jar R=${REF} O=${REF//fasta/dict}
fi

if [[ $LOCIFILE == *"bed" ]]
  then
  perl -lane 'print $F[0].":".$F[1]."-".$F[2]' $LOCIFILE > ${LOCIFILE//.bed/.intervals}
  LOCIFILE=${LOCIFILE//.bed/.intervals}
fi

echo java -jar  ${GATK}/GenomeAnalysisTK.jar \
    -T HaplotypeCaller \
    -I ${BAM} \
    -ploidy ${PLOIDY} \
    -R ${REF} \
    -o ${OUTFILE} \
    -gt_mode DISCOVERY \
    -out_mode EMIT_ALL_CONFIDENT_SITES \
    -stand_emit_conf 0 \
    -stand_call_conf 0 \
    -L ${LOCIFILE} \
    -APO ${OUTFILE/.vcf/.apo} \
    -forceActive  \
    -globalMAPQ -1 -pcrModel NONE \
    -useFilteredReadsForAnnotations \
    -mmq 0

java -jar  ${GATK}/GenomeAnalysisTK.jar \
    -T HaplotypeCaller \
    -I ${BAM} \
    -ploidy ${PLOIDY} \
    -R ${REF} \
    -o ${OUTFILE} \
    -gt_mode DISCOVERY \
    -out_mode EMIT_VARIANTS_ONLY \
    -stand_emit_conf 0 \
    -stand_call_conf 30 \
    -L ${LOCIFILE} \
    -APO ${OUTFILE/.vcf/.apo} \
    -globalMAPQ -1 -pcrModel NONE \
    -useFilteredReadsForAnnotations \
    -mmq 0 -nda
# -ERC GVCF \
# -variant_index_type LINEAR -variant_index_parameter 128000
#    -forceActive  \

#!/bin/bash

## defaults
QUERY=1
REFDIR="/seq/plasmodium/sredmond/lobSTRrefs/Pf3D7_lobstr"
RGSAMPLE="sample"
RGLIB="lib1"
RGID="lib1"
OUT="runLobSTR.out"

LOBSTR_NOISE_MODEL="/seq/plasmodium/sredmond/lobSTR/share/lobSTR/models/illumina_v3.pcrfree"


while getopts "b:r:s:l:i:t:o:h" opt; do
  case $opt in
    b) QUERY=${OPTARG} ;;
    r) REFDIR=$OPTARG ;;
    s) RGSAMPLE=$OPTARG ;;
    l) RGLIB=$OPTARG ;;
    i) RGID=$OPTARG ;;
    t) TMP=$OPTARG ;;
    o) OUT=$OPTARG ;;
   \?) echo "Invalid option: -$OPTARG" >&2 ;;
    h) echo "  runLobSTR.sh -q <query_bam> -r <ref_directory>";
       echo "              -s <readgroup_sample> -l <readgroup_lib>";
       echo "              -t <tmpdir:tmp>";;
  esac
done

TMP=$OUT

echo $QUERY
echo $OUT

REFERENCE="${REFDIR}/lobSTR_"

mkdir $TMP
cd $TMP

BAM=`basename $QUERY`
 
echo "getting sorted read-groupified file"
ln -s $QUERY $BAM
#remove secondary alignments & transfer to working folder
samtools view -b -F 256 $BAM > ${BAM/.bam/.u.bam}
#add read groups (must be sorted by pos - why??)
addReadGroups -b ${BAM/.bam/.u.bam} -o ${BAM/.bam/.rg.bam} \
    -s $RGSAMPLE -l $RGLIB -i $RGID
#resort by name for lobSTR
samtools sort -n ${BAM/.bam/.rg.bam} ${BAM/.bam/.sort}

#rm ${BAM/.bam/.rg.bam}

echo "running lobSTR align"
echo lobSTR \
   --index-prefix $REFERENCE \
   -f ${BAM/.bam/.sort.bam} --bampair \
   --rg-sample $RGSAMPLE --rg-lib $RGLIB \
   --out $OUT

lobSTR \
   --index-prefix $REFERENCE \
   -f ${BAM/.bam/.sort.bam} --bampair \
   --rg-sample $RGSAMPLE --rg-lib $RGLIB \
   --out $OUT

echo "sorting and indexing"
samtools sort ${OUT}.aligned.bam ${OUT}.sorted
samtools index ${OUT}.sorted.bam


echo "running lobSTR allelotype"

allelotype \
   --command classify \
   --bam ${OUT}.sorted.bam \
   --index-prefix $REFERENCE \
   --strinfo $REFDIR/strinfo.tab \
   --noise_model $LOBSTR_NOISE_MODEL \
   --haploid all \
   --out ${OUT}

cp ${OUT}.vcf ../

exit 0

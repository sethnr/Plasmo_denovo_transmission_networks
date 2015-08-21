#!/bin/bash

OUT=$1 
REF=$2
shift

SAMPLE_NO=`ls $OUT | wc -l`


#check if it is already there (i.e. SNPs already calculated)
ls ${OUT}_all.depth
rc=$?;
if [[ $rc == 0 ]];
then
    echo "getEdgeCoverage OUTPUT ALREADY PRESENT, NOT RUNNING";
    exit 0;
fi


if [[ $SAMPLE_NO != $VCF_NO ]]
   then
       for SAMPLE in `ls $OUT `
       do
	   if [[ -d "${OUT}/${SAMPLE}/" ]]
	       then
	       echo "concatting"
	       cat ${OUT}/${SAMPLE}/*final.fasta > ${OUT}/${SAMPLE}/${SAMPLE}.fasta
	       perl ~/bin/fasta_to_fastq.pl ${OUT}/${SAMPLE}/${SAMPLE}.fasta > ${OUT}/${SAMPLE}/${SAMPLE}.fastq
#	       bwa mem $REF ${OUT}/${SAMPLE}/${SAMPLE}.fastq | \
#		   samtools view -b - > ${OUT}/${SAMPLE}/${SAMPLE}.edges.bam
#	       samtools sort -T ${SAMPLE}.tmp -O bam -o ${OUT}/${SAMPLE}/${SAMPLE}.edges.s.bam ${OUT}/${SAMPLE}/${SAMPLE}.edges.bam

	   fi
       done
   TODEPTH=`ls ${OUT}/*/*edges.s.bam`
fi

SAMPLES=""
EDGEBAMS=""
HEADER="CHR\tPOS"
for F in `ls $OUT`
do
    if [[ -d ${OUT}/${F} ]]
	then
	echo "DIR"  ${OUT}/${F} 
	SAMPLES="${SAMPLES} ${F}"
	HEADER="${HEADER}\t${F}"
	EDGEBAMS="${EDGEBAMS} ${OUT}/${F}/${F}.edges.s.bam"
	else
	echo "FILE"  ${OUT}/${F} 
    fi
done
echo $SAMPLES

echo "depthifying"
#echo $TODEPTH
echo `ls $OUT `
samtools depth $EDGEBAMS > ${OUT}/${OUT}.edges.depth

#echo "${HEADER}" > ${OUT}/${OUT}.header.txt
#cat ${OUT}/${OUT}.header.txt ${OUT}/${OUT}.edges.depth > ${OUT}/${OUT}.edges.h.depth

perl -i -pe "print \"${HEADER}\n\" if 1..1;" ${OUT}/${OUT}.edges.depth

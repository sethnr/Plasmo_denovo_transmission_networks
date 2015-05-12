#!/bin/bash

SET=$1
LANE=$2
DATA=$3
REGION=$4

#NAME=${SET}_${STR}_${LANE}
NAME=${SET}_${LANE}

echo "SETUP..."
mkdir $NAME
cd $NAME
mkdir tmp

echo "SET " $SET
echo "LANE " $LANE
echo "DATA " $DATA
echo "REGION " $REGION

echo "LINK/INDEX BAM..."
BAMFILE=`ls ${DATA}/${LANE}/${SET}/*bam`

ln -s $BAMFILE ${NAME}.bam
samtools index ${NAME}.bam

echo "DISCOVAR"
Discovar READS=${NAME}.bam \
	 REGIONS=${REGION} \
	 OUT_HEAD=${NAME} \
	 TMP=./tmp \
	 REFERENCE=${WORK}/refs/PlasmoDB-24_Pfalciparum3D7_Genome.fasta

echo "MAKE OUTPUTS"
dot -Tpng -o ${NAME}.final.png ${NAME}.final.dot
perl -i -pe "s/${SET}$/${NAME}/gi" ${NAME}.final.variant.filtered.vcf
bgzip ${NAME}.final.variant.filtered.vcf
tabix ${NAME}.final.variant.filtered.vcf.gz

echo "CLEANUP"
rm -r tmp



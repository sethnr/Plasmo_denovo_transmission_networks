#!/bin/bash

SET=$1
LANE=$2
STR=$3
REGION=$4

#NAME=${SET}_${STR}_${LANE}
NAME=${SET}_${LANE}

mkdir $NAME
cd $NAME
mkdir tmp

FILE=`ls ${DISCODATA}/${LANE}/${SET}/*bam`

ln -s $FILE ${NAME}.bam

Discovar READS=${NAME}.bam \
	 REGIONS=${REGION} \
	 OUT_HEAD=${NAME} \
	 TMP=./tmp \
	 REFERENCE=${WORK}/refs/PlasmoDB-24_Pfalciparum3D7_Genome.fasta

echo "DISCOVAR COMPLETE"
rm -r tmp

dot -Tpng -o ${NAME}.final.png ${NAME}.final.dot

perl -i -pe "s/${SET}$/${NAME}/gi" ${NAME}.final.variant.filtered.vcf
bgzip ${NAME}.final.variant.filtered.vcf
tabix ${NAME}.final.variant.filtered.vcf.gz



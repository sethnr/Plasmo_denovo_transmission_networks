#!/bin/bash

SET=$1
LANE=$2
STR=$3
REGION=$4

NAME=${SET}_${STR}_${LANE}

mkdir $NAME
cd $NAME
mkdir tmp

FILE=`ls ${DISCODATA}/${LANE}/${SET}/*bam`

Discovar READS=${FILE} \
	 REGIONS=${REGION} \
	 OUT_HEAD=${NAME} \
	 TMP=./tmp \
	 REFERENCE=${WORK}/refs/PlasmoDB-24_Pfalciparum3D7_Genome.fasta

rm -r tmp

dot -Tpng -o ${NAME}.final.png ${NAME}.final.dot

#!/usr/bin/bash

Use Discovar

set SET=$1
set LANE=$2
set STR=$3
set REGION=$4

set NAME=${SET}_${STR}_${LANE}

mkdir $NAME
cd $NAME
mkdir tmp

set FILE=`ls ${DISCODATA}/${LANE}/${SET}/*bam`

Discovar READS=${FILE} \
	 REGIONS=${REGION} \
	 OUT_HEAD=${NAME} \
	 TMP=./tmp
	 REFERENCE=${WORK}/refs/PlasmoDB-24_Pfalciparum3D7_Genome.fasta

rm -r tmp

dot -Tpng -o ${NAME}.final.png ${NAME}.final.dot

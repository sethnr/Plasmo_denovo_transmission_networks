#!/bin/bash

GENOME=$1
NAME=`basename $GENOME`
NAME=${NAME/.fasta/}

GENOME=`readlink -f $GENOME`
mkdir ${NAME}_split
cd ${NAME}_split
pyfasta split --header="%(seqid)s.fasta"  $GENOME

for FASTA in *fasta 
do
  FASTA=${FASTA/.fasta/}
  $DISCO1/STRcallers/mreps/mreps -fasta ${FASTA}.fasta \
     | perl -lane "print join(\"\t\",${FASTA},@F) unless 1..15" \
     | grep -e '->' 
done  \
> ../STRs_mreps_${NAME}.R0.txt

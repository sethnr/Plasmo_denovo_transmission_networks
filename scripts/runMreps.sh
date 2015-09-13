#!/bin/bash

GENOME=$1
R=$2

if [[ -z "$R" ]]; then 
R=0
fi

NAME=`basename $GENOME`
NAME=${NAME/.fasta/}

GENOME=`readlink -f $GENOME`
mkdir ${NAME}_split
cd ${NAME}_split
pyfasta split --header="%(seqid)s.fasta"  $GENOME

for FASTA in *fasta 
do
  FASTA=${FASTA/.fasta/}
#  echo $DISCO1/STRcallers/mreps/mreps -res $R -fasta ${FASTA}.fasta
  $DISCO1/STRcallers/mreps/mreps -res $R -fasta ${FASTA}.fasta \
     | perl -lane "print join(\"\t\",${FASTA},@F) unless 1..15" \
     | grep -e '->' 
done  \
> ../STRs_mreps_${NAME}.R${R}.txt

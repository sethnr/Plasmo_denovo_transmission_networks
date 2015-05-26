#!/bin/bash


## defaults
REF=$WORK/refs/Pf3D7_v3.fasta
OUTDIR="alignViewOut"

while getopts "f:p:o:" opt; do
  case $opt in
      r) REF=${OPTARG};;
      p) REGION=${OPTARG};;
      o) OUTDIR=${OPTARG};;
   \?) echo "Invalid option: -$OPTARG" >&2 ;;
  esac
done
shift $((OPTIND-1))

OPTIND=1

echo "aligning contigs with bwmem"

if [ ! -d ${OUTDIR} ];
then
  mkdir ${OUTDIR}
fi

cd $OUTDIR

for FASTA in "$@"
do
	perl  ~/bin/fasta_to_fastq.pl ../$FASTA > ${FASTA//.fasta/.fastq}
	echo bwa mem  -SP $REF ${FASTA//.fasta/.fastq}
	bwa mem  -SP $REF ${FASTA//.fasta/.fastq} \
	 > ${FASTA//.fasta/.sam}
	samindex ${FASTA//.fasta/.sam}
done

echo "new" > viewRegion.${REGION}.igv
echo "genome "${REF} >> viewRegion.${REGION}.igv
echo "goto "${REGION} >> viewRegion.${REGION}.igv

for FASTA in "$@"
do
 echo "load " ${FASTA//.fasta/.bam} >> viewRegion.${REGION}.igv
done

echo "collapse" >> viewRegion.${REGION}.igv
echo "snapshot "${OUTDIR}".png" >> viewRegion.${REGION}.igv

#use igv-2.3

igv.sh -b viewRegion.${REGION}.igv


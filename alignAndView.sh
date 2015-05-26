#!/bin/bash


## defaults
REF=$WORK/refs/Pf3D7_v3.fasta

while getopts "f:p:" opt; do
  case $opt in
      r) REF=${OPTARG}
	 shift $((OPTIND-1));;
      p) REGION=${OPTARG}
	 shift $((OPTIND-1));;
   \?) echo "Invalid option: -$OPTARG" >&2 ;;
  esac
done
OPTIND=1

echo "aligning contigs with bwmem"

for FASTA in "$@"
do
	perl  ~/bin/fasta_to_fastq.pl $FASTA > ${FASTA//.fasta/.fastq}
	bwa mem  -SP  $REFERENCE \
	 ${FASTA//.fasta/.fastq} \
	 > ${FASTA//.fasta/.sam}
	samindex ${FASTA//.fasta/.sam}
done

cat "new" > viewRegion.${REGION}.igv
cat "genome "${REF} >> viewRegion.${REGION}.igv
cat "goto "${REGION} >> viewRegion.${REGION}.igv

for FASTA in "$@"
do
 echo "load " ${FASTA//.fasta/.bam} >> viewRegion.${REGION}.igv
done

cat "snapshot" >> viewRegion.${REGION}.igv

use igv-2.3

igv -b viewRegion.${REGION}.igv

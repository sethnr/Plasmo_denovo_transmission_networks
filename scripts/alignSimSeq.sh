#!/bin/bash


## defaults
QUERY=1
REFERENCE=2
READFRAGSIZE=450
READPAIRS=60000 #~100x for Pf3D7 genome
READLENGTH=250
OUT="alignSimSeq.out"

while getopts "q:r:s:n:l:o:h" opt; do
  case $opt in
    q) QUERY=${OPTARG} ;;
    r) REFERENCE=$OPTARG ;;
    s) READFRAGSIZE=$OPTARG ;;
    n) READPAIRS=$OPTARG ;;
    l) READLENGTH=$OPTARG;;
    t) TMP=$OPTARG ;;
    o) OUT=$OPTARG ;;
   \?) echo "Invalid option: -$OPTARG" >&2 ;;
    h) echo "  alignSimSeq -q <query_fasta> -r <ref_fasta>";
       echo "              -s <fragment_size:450> -l <read_length:250>";
       echo "              -n <no_read_pairs:60,000> -t <tmpdir:tmp>";;
  esac
done

QUERYOUT=${QUERY/%.gz/}
QUERYOUT=${QUERYOUT/%.fasta/.fastq}
QUERYOUT=${QUERYOUT/%.fa/.fastq}
QOUT1=${QUERYOUT/%.fastq/.r1.fastq}
QOUT2=${QUERYOUT/%.fastq/.r2.fastq}


echo wgsim -r 0 \
-N $READPAIRS -1${READLENGTH} -2${READLENGTH} \
$QUERY  $QOUT1 $QOUT2
wgsim -r 0 \
-N $READPAIRS -1${READLENGTH} -2${READLENGTH} \
$QUERY  $QOUT1 $QOUT2

if [[ ! -f ${QOUT1} &&  -f ${QOUT2} ]];
then
    echo "wgsim not run?"
    exit 1
fi


echo bwa mem -M -t 5 $REFERENCE  $QOUT1  $QOUT2 \| samtools view -bS - \> ${OUT}.bam
bwa mem -M -t 5 $REFERENCE  $QOUT1  $QOUT2 | samtools view -bS - > ${OUT}.bam

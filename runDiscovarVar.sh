#!/bin/bash

#SET=$1
#LANE=$2
#DATA=$3
#REGION=$4

DATA="./"
$OPTIND=0;
while getopts ":a" opt; do
  case $opt in
    s) SET=$OPTARG ;;
    l) LANE=$OPTARG ;;
    d) DATA=$OPTARG ;;
    r) REGION=$OPTARG ;;
    m) MEMORY=$OPTARG ;;
    \?) echo "Invalid option: -$OPTARG" >&2 ;;
  esac
done

SET=${@:$OPTIND:1}
LANE=${@:$OPTIND+1:1}
REGION=${@:$OPTIND+2:1}



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

ls ${NAME}.final.variant.filtered.vcf
rc=$?;
if [[ $rc == 0 ]];
then
    echo "DISCOVAR OUTPUT ALREADY PRESENT, NOT RUNNING";
    exit 0;
fi


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

rc=$?;
if [[ $rc != 0 ]];
then
    echo "DISCOVAR ERROR";
    exit $rc;
fi

ls ${NAME}.final.variant.filtered.vcf
rc=$?;
if [[ $rc != 0 ]];
then
    echo "DISCOVAR FINISHED WITHOUT OUTPUT";
    exit $rc;
fi


echo "MAKE OUTPUTS"
dot -Tpng -o ${NAME}.final.png ${NAME}.final.dot
perl -i -pe "s/${SET}$/${NAME}/gi" ${NAME}.final.variant.filtered.vcf
bgzip ${NAME}.final.variant.filtered.vcf
tabix ${NAME}.final.variant.filtered.vcf.gz

if [[ $rc != 0 ]];
then
    echo "TABIX/BGZIP ERROR";
    exit $rc;
fi

echo "CLEANUP"
rm -r tmp



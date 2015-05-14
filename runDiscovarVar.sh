#!/bin/bash

#SET=$1
#LANE=$2
#DATA=$3
#REGION=$4

#VARIABLE=<< EOL
#./runDiscovarVar.sh [options] library  lane  region,region...
#options:
#     -d /path/to/data/directory
#     -f ./reference.fasta
#     -m max_mem (GB) [2]
#     -n nodes [1]
#EOL

#DATA="./"
#OPTIND=1

## defaults
NODES=1
MEMORY=2
SPLIT=0

while getopts "d:f:m:n:s" opt; do
  case $opt in
    d) DATA=${OPTARG} ;;
    f) REFERENCE=$OPTARG ;;
    m) MEMORY=$OPTARG ;;
    n) NODES=$OPTARG ;;
    s) SPLIT=1;;
#    s) SET=$OPTARG ;;
#    l) LANE=$OPTARG ;;
#    r) REGION=$OPTARG ;;
   \?) echo "Invalid option: -$OPTARG" >&2 ;;
  esac
done

echo "OPTIND - " $OPTIND

SET=${@:$OPTIND:1}
LANE=${@:$OPTIND+1:1}
REGION=${@:$OPTIND+2:1}


#only get one bam for all regions
#put all regions in same file
REFNAME=${SET}_${LANE}
#IF NOT SPLITTING REGIONS, THEN NAME = SET+NAME
NAME=${SET}_${LANE}
#OTHERWISE ADD REGION TO END
if [[ $SPLIT==1 ]]; then NAME=${NAME}_$REGION; fi

#if [[ $SPLIT==1 ]]; then REFNAME=${NAME}_$REGION; fi
#NB, inefficient - for the moment, make new bam index for each bam, 
# even though most are the same file - avoids race condition across farm
# eventually should make NAME_LANE.bam.bai


echo "SETUP..."
mkdir $REFNAME
cd $REFNAME
mkdir tmp


echo "SET " $SET
echo "LANE " $LANE
echo "DATA " $DATA
echo "REGION " $REGION
echo "MEMORY " $MEMORY
echo "NODES " $NODES

# DO SOME STUFF:

ls ${NAME}.final.variant.filtered.vcf
rc=$?;
if [[ $rc == 0 ]];
then
    echo "DISCOVAR OUTPUT ALREADY PRESENT, NOT RUNNING";
    exit 0;
fi


echo "LINK/INDEX BAM..."
BAMFILE=`ls ${DATA}/${LANE}/${SET}/*bam`

#check if bam index file already present
ls ${NAME}.bam.bai
rc=$?;
if [[ $rc != 0 ]];
then
    ln -s $BAMFILE ${NAME}.bam
    samtools index ${NAME}.bam
fi

echo "DISCOVAR"
Discovar READS=${NAME}.bam \
	 REGIONS=${REGION} \
	 OUT_HEAD=${NAME} \
	 TMP=./tmp \
	 REFERENCE=${WORK}/refs/PlasmoDB-24_Pfalciparum3D7_Genome.fasta \
         NUM_THREADS=${NODES} \
	 MAX_MEMORY_GB=${MEMORY}
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



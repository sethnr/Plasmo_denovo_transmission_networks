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
REFERENCE=${WORK}/refs/PlasmoDB-24_Pfalciparum3D7_Genome.fasta \

while getopts "d:f:m:n:s:l:B:N:" opt; do
  case $opt in
    d) DATA=${OPTARG} ;;
    f) REFERENCE=$OPTARG ;;
    m) MEMORY=$OPTARG ;;
    n) NODES=$OPTARG ;;
    s) SPLIT=1;;
    B) BAMFILE=$OPTARG;;
    N) NAME=$OPTARG;;
#    s) SET=$OPTARG ;;
    l) LANE=$OPTARG ;;
#    r) REGION=$OPTARG ;;
   \?) echo "Invalid option: -$OPTARG" >&2 ;;
  esac
done

echo "OPTIND - " $OPTIND

if [ -z "$BAMFILE" ]; then
    SET=${@:$OPTIND:1}
    LANE=${@:$OPTIND+1:1}
    REGION=${@:$OPTIND+2:1}
#only get one bam for all regions
#put all regions in same file
    REFNAME=${SET}_${LANE}
#IF NOT SPLITTING REGIONS, THEN NAME = SET+NAME
    NAME=${SET}_${LANE}
#OTHERWISE ADD REGION TO END
    BAMFILE=`ls ${DATA}/${LANE}/${SET}/*bam`

elif [ -z "$NAME" ]; then
     echo "NAME [-N] must be set if BAMFILE [-B] is set"
     exit 1
#elif [ -z "$LANE" ]; then 
else:
     echo "BAMFILE set: "${BAMFILE}
     echo "NAME set: "${NAME}
     REFNAME=$NAME
     SET=$NAME
     REGION=${@:$OPTIND:1}
#else
#     echo "BAMFILE set: "${BAMFILE}
#     echo "NAME set: "${NAME}_${LANE}
#     REFNAME=${NAME}_${LANE}
#     SET=$NAME
#     REGION=${@:$OPTIND:1}
fi


if [[ $SPLIT==1 ]]; then NAME=${NAME}_$REGION; fi

#if [[ $SPLIT==1 ]]; then REFNAME=${NAME}_$REGION; fi
#NB, inefficient - for the moment, make new bam index for each bam, 
# even though most are the same file - avoids race condition across farm
# eventually should make NAME_LANE.bam.bai


echo "SETUP..."
mkdir $REFNAME
cd $REFNAME
#mkdir tmp_${NAME}
echo "SET " $SET
echo "LANE " $LANE
echo "DATA " $DATA
echo "REGION " $REGION
echo "MEMORY " $MEMORY
echo "NODES " $NODES

# DO SOME STUFF:

#ls ${NAME}.final.variant.filtered.vcf
#rc=$?;
#if [[ $rc == 0 ]];

#if [[ -f ${NAME}.final.variant.filtered.vcf ]]
#then
#    echo "DISCOVAR OUTPUT ALREADY PRESENT, NOT RUNNING";
#    exit 0;
if [[ -f ${NAME}.final.variant.filtered.vcf.gz ]]
then
    echo "DISCOVAR OUTPUT ALREADY PRESENT, NOT RUNNING";
    exit 0;
fi

#MAKE SUB FOLDER (FULL RUN IS MAKING TOO MANY FILES IN ONE FOLDER)
mkdir ${NAME}_${LANE}
cd ${NAME}_${LANE}
mkdir tmp_${NAME}

echo "WORKING DIR"
pwd

if [[ ! -f ${NAME}.final.variant.filtered.vcf ]];
then
    echo "LINK/INDEX BAM..."

#check if bam index file already present
##    ls ${NAME}.bam.bai
##    rc=$?;
##    if [[ $rc != 0 ]];
##	then
##	ln -s $BAMFILE ${NAME}.bam
##	echo samtools index ${NAME}.bam
##	samtools index ${NAME}.bam
##    fi
    #remove old bams - some are wrong
    if [[ -e ${NAME}.bam ]]; then
	echo "REMOVING OLD BAM"
	rm ${NAME}.bam;
	rm ${NAME}.bam.bai
    fi

    ln -s $BAMFILE ${NAME}.bam
    echo samtools index ${NAME}.bam
    samtools index ${NAME}.bam
    
#mkdir tmp_${NAME}

    echo "DISCOVAR"
    echo Discovar READS=${NAME}.bam \
	REGIONS=${REGION} \
	OUT_HEAD=${NAME} \
	TMP=./tmp_${NAME} \
	REFERENCE=${REFERENCE} \
	NUM_THREADS=${NODES} \
	 MAX_MEMORY_GB=${MEMORY}
    Discovar READS=${NAME}.bam \
	REGIONS=${REGION} \
	OUT_HEAD=${NAME} \
	TMP=./tmp_${NAME} \
	REFERENCE=${REFERENCE} \
	NUM_THREADS=${NODES} \
	 MAX_MEMORY_GB=${MEMORY}
    rc=$?;
    if [[ $rc == 1 ]];
	then
	echo "DISCOVAR ERROR";
	#catch instances where finished with no output
	#why doesn't disco output 0?
	if [[ -f ${NAME}.final.fasta ]]
	    then
	    ls -l ${NAME}.final*
	    echo "DISCOVAR FINISHED WITH ERROR CODE",$rc;
	    exit 0
	fi
	exit $rc;
    elif [[ $rc != 0 ]];
	then
	echo "DISCOVAR ERROR CODE",$rc;
	exit $rc;
    fi

    #if error code is zero:
    if [[ ! -f ${NAME}.final.variant.filtered.vcf ]]
	then
	echo "DISCOVAR FINISHED, BUT WITHOUT OUTPUT";
	exit 0
    fi
fi


echo "MAKE OUTPUTS"
dot -Tpng -o ${NAME}.final.png ${NAME}.final.dot
# perl -i -pe "s/${SET}$/${NAME}/gi" ${NAME}.final.variant.filtered.vcf

#if using set/lane/etc, change back to refname before concatenating
if [[ -n "$SET" && -n "$LANE" ]]; then
    perl -i -pe "s/${SET}$/${REFNAME}/gi" ${NAME}.final.variant.filtered.vcf
fi 

#UGLY HACK (discovarRegion is left blank)
perl -i -ne 'print $_ unless $_ =~ m/DiscovarRegion/gi' ${NAME}.final.variant.filtered.vcf
#UGLY HACK sample name is 'unknown' in fakeNGS samples
perl -i -pe "s/unknown(?>$)/${SET}/" ${NAME}.final.variant.filtered.vcf

#cleanup:
rm ${NAME}.{0..9}*

#move back to parent folder
cp ${NAME}.final.variant.filtered.vcf ../${NAME}.final.variant.filtered.vcf
cd ../

bgzip ${NAME}.final.variant.filtered.vcf
tabix ${NAME}.final.variant.filtered.vcf.gz

$tabrc=$?
if [[ $tabrc != 0 ]];
then
    echo "TABIX/BGZIP ERROR";
    exit $rc;
fi

echo "CLEANUP"
# rm -r tmp_${NAME}

exit 0

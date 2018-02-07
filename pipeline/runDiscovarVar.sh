#!/bin/bash

source /broad/software/scripts/useuse

reuse -q Samtools
reuse -q BWA
reuse -q .vcftools-0.1.14
reuse -q Tabix
reuse -q Python-2.7
reuse -q .numpy-1.9.1-python-2.7.1-sqlite3-rtrees

reuse -q GCC-4.9
reuse -q Discovar


## defaults
NODES=1
MEMORY=2
SPLIT=0
REFERENCE=${WORK}/refs/PlasmoDB-24_Pfalciparum3D7_Genome.fasta \


	 
while getopts "f:m:n:B:o:r:" opt; do
  case $opt in
    f) REFERENCE=$OPTARG ;;
    m) MEMORY=$OPTARG ;;
    n) NODES=$OPTARG ;;
    B) BAMLIST+=($OPTARG);;
    o) NAME=$OPTARG;;
    r) REGION=$OPTARG ;;
   \?) echo "Invalid option: -$OPTARG" >&2 ;;
  esac
done


echo "SETUP..."
mkdir $NAME

echo "BAM " $BAM
echo "REGION " $REGION
echo "MEMORY " $MEMORY
echo "NODES " $NODES


for BAM in "${BAMLIST[@]}"; do
#    ln -s $BAM
#    ln -s ${BAM}.idx    
    #    BAMSTR=${BAMSTR}"  READS="`basename $BAM`
    ln -s `readlink -e $BAM` ${NAME}/`basename $BAM`
    ln -s `readlink -e ${BAM}.bai` ${NAME}/`basename ${BAM}.bai`
    BAMLIST2+=(`basename $BAM`)
done

cd $NAME
mkdir tmp_$NAME


#BAMSTR=shift  "${BAMLIST[@]}"
BAMSTR=$(IFS=, ; echo "${BAMLIST2[*]}")



# DO SOME STUFF:
if [[ -f ${NAME}.final.variant.filtered.vcf.gz ]]
then
    echo "DISCOVAR OUTPUT ALREADY PRESENT, NOT RUNNING";
    exit 0;
fi


echo "WORKING DIR"
pwd


echo "DISCOVAR"
echo Discovar READS=$BAMSTR \
	REGIONS=${REGION} \
	OUT_HEAD=${NAME} \
	TMP=./tmp_${NAME} \
	REFERENCE=${REFERENCE} \
	NUM_THREADS=${NODES} \
	 MAX_MEMORY_GB=${MEMORY}
Discovar READS=$BAMSTR \
	REGIONS=${REGION} \
	OUT_HEAD=${NAME} \
	TMP=./tmp_${NAME} \
	REFERENCE=${REFERENCE} \
	NUM_THREADS=${NODES} \
	 MAX_MEMORY_GB=${MEMORY}
rc=$?;
echo "exit" $rc
exit 0 

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

tabrc=$?
if [[ $tabrc -ne 0 ]];
then
    echo "TABIX/BGZIP ERROR";
    exit $rc;
fi

echo "CLEANUP"
# rm -r tmp_${NAME}

exit 0

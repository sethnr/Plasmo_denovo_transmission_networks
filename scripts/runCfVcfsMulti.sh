#!/bin/bash

##join function (use to pass args to python script)
##function join { local IFS="$1"; shift; echo "$*"; }

DEPTH=$DISCO1/disco_halfmeg_fakeNGS/fakev3D7.depth.gz

while getopts "b:r:o:c:s:e:" opt; do
  case $opt in
    b) BED=${OPTARG} ;;
    r) REGION=$OPTARG ;;
    o) OUTFILE=$OPTARG ;;
    c) CHROM=$OPTARG ;;
    s) START=$OPTARG ;;
    e) END=$OPTARG ;;
    d) DEPTH=$OPTARG ;;
   \?) echo "Invalid option: -$OPTARG" >&2 ;;
    h) echo "  runCfVcfsMulti.sh -r <region> -b <bed file> (for filtering matched regions)";
       echo "              -o <outfile> vcffile1 vcffile2 vcffileN";;
  esac
done
shift $(expr $OPTIND - 1 )
echo $@


echo "setup"

#get abspaths
BED=`readlink -f  $BED`

echo "bed" $BED
echo "region" $REGION
echo "outfile" $OUTFILE

mkdir tmp_${OUTFILE}
TMP=tmp_${OUTFILE}

#make local links
#ln -s $BED ${TMP}/`basename $BED`
#BED=${TMP}/`basename $BED`
ln -s $BED ${TMP}/`basename $BED`
BED=${TMP}/`basename $BED`

VCFS=()
BASEVCFS=()
echo "parsing matching regions"
for VCF in $@
do
  echo $VCF
  VCFARR=(${VCF/:::/ })
  echo ${VCFARR[0]}
  echo ${VCFARR[1]}
#  BASEVCF=`basename ${VCF} `
 
  if [[ ${#VCFARR[@]} > 1 ]]
    then
      echo "SET"
      #split into name/vcf 
      VCF=${VCFARR[1]}
      NAME=${VCFARR[0]}

      VCF=`readlink -f $VCF`
      BASEVCF=${NAME}.vcf.gz

      #change to do everything in tmp folder
      VCFS+=(${NAME}:::${BASEVCF})
#      VCFS+=(${NAME}:::${TMP}/${BASEVCF})
      BASEVCFS+=($BASEVCF)
    else
      echo "NOT SET"
      VCF=`readlink -f $VCF`
      BASEVCF=`basename ${VCF} `
      NAME=$BASEVCF
      VCFS+=(${BASEVCF})
#      VCFS+=(${TMP}/${BASEVCF})
      BASEVCFS+=($BASEVCF)
  fi
  echo $BASEVCF
#  VCF=`readlink -f $VCF`
  if [[ -n "$CHROM" ]]
    then 
      vcftools --gzvcf $VCF --recode --out ${TMP}/${BASEVCF/.vcf.gz/} --bed $BED \
        --chr $CHROM --from-bp $START --to-bp $END
#      vcftools --gzvcf $VCF --recode --out ${BASEVCF/.vcf.gz/} --bed $BED \
#        --chr $CHROM --from-bp $START --to-bp $END
    else
      vcftools --gzvcf $VCF --recode --out ${TMP}/${BASEVCF/.vcf.gz/} --bed $BED
#      vcftools --gzvcf $VCF --recode --out ${BASEVCF/.vcf.gz/} --bed $BED
  fi
done

echo "cd into tmp for analysis"
cd $TMP
echo "bgzipping everything"
for BASEVCF in ${BASEVCFS[@]}
do
#  mv ${TMP}/${BASEVCF/.vcf.gz/.recode.vcf} ${TMP}/${BASEVCF/.vcf.gz/.vcf}
#  bgzip ${TMP}/${BASEVCF/.vcf.gz/.vcf}
#  tabix -pvcf ${TMP}/${BASEVCF}
  mv ${BASEVCF/.vcf.gz/.recode.vcf} ${BASEVCF/.vcf.gz/.vcf}
  bgzip ${BASEVCF/.vcf.gz/.vcf}
  tabix -pvcf ${BASEVCF}
# VCFS+=($BASEVCF)
done


#echo ${VCFS[@]}
#THIS CAN'T BE THE BEST WAY TO DO THIS:
VCFLINE=`for i in "${VCFS[@]}" ; do echo "-v ${i}" ; done`
#echo $VCFLINE
echo "comparing vcfs"
echo python $DISCO1/scripts/cfs/cfVCFsMulti.py \
  -o ${OUTFILE}.vcf \
  $VCFLINE
#  -r $REGION \

python $DISCO1/scripts/cfs/cfVCFsMulti.py \
  -o ${OUTFILE}.vcf \
  $VCFLINE
#  -r $REGION \

echo "adding STR information"
echo python $DISCO1/scripts/cfs/addSTRsToVCF.py -v ${OUTFILE}.vcf -s $DISCO1/analyses/STR_calling/STRs_mreps_Pf3D7_v3.R0.txt
python $DISCO1/scripts/cfs/addSTRsToVCF.py -v ${OUTFILE}.vcf -s $DISCO1/analyses/STR_calling/STRs_mreps_Pf3D7_v3.R0.txt


echo "calculating match statistics"
echo python $DISCO1/scripts/cfs/getVcfMatchStats.py -v ${OUTFILE}.STRs.vcf \> ${OUTFILE}.txt
python $DISCO1/scripts/cfs/getVcfMatchStats.py -v ${OUTFILE}.STRs.vcf > ${OUTFILE}.txt

echo "calculating depth/match statistics"
echo python $DISCO1/scripts/cfs/getDepthHetMulti.py -v ${OUTFILE}.STRs.vcf \
    -d $DEPTH \
    -b 5 -r ${CHROM}:${START}-${END} \> ${OUTFILE}.DEPTH.txt
python $DISCO1/scripts/cfs/getDepthHetMulti.py -v ${OUTFILE}.STRs.vcf \
    -d $DEPTH \
    -b 5 -r ${CHROM}:${START}-${END} > ${OUTFILE}.DEPTH.txt

cp ${OUTFILE}.txt ../${OUTFILE}.txt
cp ${OUTFILE}.DEPTH.txt ../${OUTFILE}.DEPTH.txt
cp ${OUTFILE}.STRs.vcf.gz ../${OUTFILE}.vcf.gz

exit 0

#!/bin/bash

COMPARISON=$1
SAMPLE=$2
OUTFILE=$3
MINOR=$4

#minor = '-M' - passed directly to relevant scripts

SAMPLE_ABS=`readlink -f ${SAMPLE}`
COMPARISON_ABS=`readlink -f ${COMPARISON}`
SAMPLE=`basename $SAMPLE`
COMPARISON=`basename $COMPARISON`

if [[ -z $OUTFILE ]] 
then
    OUTFILE=$SAMPLE
fi

mkdir tmp_${OUTFILE}_cfDD2
cd tmp_${OUTFILE}_cfDD2

ln -s $SAMPLE_ABS
ln -s $COMPARISON_ABS

REF="${WORK}/refs/Pf3D7_v3.fasta"
GATK=/humgen/gsa-hpprojects/GATK/bin/GenomeAnalysisTK-3.4-0-g7e26428


echo "CHECK SAMPLE CONSISTENCY DD2ALL PASS"
echo python $DISCO1/scripts/cfs_dd2/checkSampleConsistencyVCF.py $MINOR \
    -v ${SAMPLE} -p FDK2D4 \
    -s SM-7LV8O -s SM-7LV8P -s SM-7LV8Q  -s SM-7LV8R
python $DISCO1/scripts/cfs_dd2/checkSampleConsistencyVCF.py $MINOR \
    -v ${SAMPLE} -p FDK2D4 \
    -s SM-7LV8O -s SM-7LV8P -s SM-7LV8Q  -s SM-7LV8R
if [[ $? != 0 ]]; then exit 1 ; fi

echo "DD2FDK V REF (PFX) DD2 PASS ONLY" 
python $DISCO1/scripts/cfs_dd2/checkRealVFakeDD2.py \
    -v1 ${SAMPLE/.vcf/.FDK2D4.vcf} -v2 $COMPARISON $MINOR \
    -s SM-7LV8O -s SM-7LV8P -s SM-7LV8Q  \
    -p FDK \
    > ${SAMPLE/.vcf/.FDK2D4.FDKCONC.txt}
if [[ $? != 0 ]]; then exit 1 ; fi

echo "DD22D4 V REF (PFX) DD2 PASS ONLY" 
python $DISCO1/scripts/cfs_dd2/checkRealVFakeDD2.py \
    -v1 ${SAMPLE/.vcf/.FDK2D4.FDKCONC.vcf} -v2 $COMPARISON $MINOR \
    -s SM-7LV8R \
    -p 2D4 \
    > ${SAMPLE/.vcf/.FDK2D4.FDKCONC.2D4CONC.txt}
if [[ $? != 0 ]]; then exit 1 ; fi

echo "CHECK SAMPLE CONSISTENCY 3D7DD2 PASS"
echo python $DISCO1/scripts/cfs_dd2/checkSampleConsistencyVCF.py $MINOR \
    -v ${SAMPLE/.vcf/.FDK2D4.FDKCONC.2D4CONC.vcf} -p 3D72D4 \
    -s SM-7LV8S -s SM-7LV8T -s SM-7LV8U  -s SM-7LV8R
python $DISCO1/scripts/cfs_dd2/checkSampleConsistencyVCF.py $MINOR \
    -v ${SAMPLE/.vcf/.FDK2D4.FDKCONC.2D4CONC.vcf} -p 3D72D4 \
    -s SM-7LV8S -s SM-7LV8T -s SM-7LV8U  -s SM-7LV8R
if [[ $? != 0 ]]; then exit 1 ; fi

echo "CHECK SAMPLE CONSISTENCY 3D7DD2 PASS"
echo python $DISCO1/scripts/cfs_dd2/checkSampleConsistencyVCF.py $MINOR \
    -v ${SAMPLE/.vcf/.FDK2D4.FDKCONC.2D4CONC.3D72D4.vcf} -p 3D7FDK \
    -s SM-7LV8S -s SM-7LV8T -s SM-7LV8U  -s SM-7LV8R
python $DISCO1/scripts/cfs_dd2/checkSampleConsistencyVCF.py $MINOR \
    -v ${SAMPLE/.vcf/.FDK2D4.FDKCONC.2D4CONC.3D72D4.vcf} -p 3D7FDK \
    -s SM-7LV8S -s SM-7LV8T -s SM-7LV8U  -s SM-7LV8R
if [[ $? != 0 ]]; then exit 1 ; fi


#echo "REAL V FAKE DD2 ALL" 
#python $DISCO1/scripts/cfs_dd2/checkRealVFakeDD2.py \
#    -v1 ${SAMPLEPOST} -v2 $COMPARISON $MINOR \
#    -s SM-7LV8O -s SM-7LV8P -s SM-7LV8Q  -s SM-7LV8R \
#    > ${SAMPLE/.vcf/.DD2CONC.txt}


#cp ${SAMPLE/.vcf/.PASS.DD2CONC.vcf} ../${OUTFILE}.PASS.DD2CONC.vcf
cp ${SAMPLE/.vcf/.FDK2D4.FDKCONC.2D4CONC.vcf} ../${OUTFILE}.PASS.DD2CONC.PFX.vcf


# python $DISCO1/scripts/cfs/addSTRsToVCF.py -s ../STRs_mreps_Pf3D7_v3.R0.txt -v ../${OUTFILE}.PASS.DD2CONC.vcf

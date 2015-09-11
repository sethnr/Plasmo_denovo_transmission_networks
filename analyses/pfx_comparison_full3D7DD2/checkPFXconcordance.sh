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

echo "FILTER HYPERVARIABLES"
python $DISCO1/scripts/cfs_dd2/filterHypervariableRecordsVCF.py -v ${SAMPLE}

echo "CHECK LANE CONCORDANCE DD2"
echo python $DISCO1/scripts/cfs_dd2/checkLaneConcordance.py \
    -v ${SAMPLE/.vcf/.HYPF.vcf} -p DD2 \
    -s SM-7LV8O -s SM-7LV8P -s SM-7LV8Q -s SM-7LV8R
python $DISCO1/scripts/cfs_dd2/checkLaneConcordance.py \
    -v ${SAMPLE/.vcf/.HYPF.vcf} -p DD2 \
    -s SM-7LV8O -s SM-7LV8P -s SM-7LV8Q -s SM-7LV8R
if [[ $? != 0 ]]; then exit 1 ; fi

echo "CHECK LANE CONCORDANCE 3D7"
echo python $DISCO1/scripts/cfs_dd2/checkLaneConcordance.py \
    -v ${SAMPLE/.vcf/.HYPF.DD2LCHK.vcf} -p 3D7 \
    -s SM-7LV8S -s SM-7LV8T -s SM-7LV8U -s SM-7LV8V
python $DISCO1/scripts/cfs_dd2/checkLaneConcordance.py \
    -v ${SAMPLE/.vcf/.HYPF.DD2LCHK.vcf} -p 3D7 \
    -s SM-7LV8S -s SM-7LV8T -s SM-7LV8U -s SM-7LV8V
if [[ $? != 0 ]]; then exit 1 ; fi

echo "CHECK SAMPLE CONSISTENCY DD2FDK"
echo python $DISCO1/scripts/cfs_dd2/checkSampleConsistencyVCF.py \
    -v ${SAMPLE/.vcf/.HYPF.DD2LCHK.3D7LCHK.vcf} -p DD2FDK  $MINOR \
    -s SM-7LV8O -s SM-7LV8P -s SM-7LV8Q
python $DISCO1/scripts/cfs_dd2/checkSampleConsistencyVCF.py \
    -v ${SAMPLE/.vcf/.HYPF.DD2LCHK.3D7LCHK.vcf} -p DD2FDK  $MINOR \
    -s SM-7LV8O -s SM-7LV8P -s SM-7LV8Q
if [[ $? != 0 ]]; then exit 1 ; fi


SAMPLEPOST=${SAMPLE/.vcf/.HYPF.DD2LCHK.3D7LCHK.DD2FDK.vcf}

#count filters
echo "FILTERING:"

perl $DISCO1/scripts/vcf/countFilters.pl ${SAMPLEPOST} >  ${SAMPLE/.vcf/.failtab.comb.txt}

#remove filtered and non-varying loci
#remove non-varying alleles (hypervariables?)
java -jar $GATK/GenomeAnalysisTK.jar \
   -T SelectVariants \
   -R $REF \
   -V ${SAMPLEPOST} \
   -o ${SAMPLE/.vcf/.PASS.vcf} \
   -env \
   -ef \
   -trimAlternates

#echo "CHECK SAMPLE CONSISTENCY DD2ALL"
#echo python $DISCO1/scripts/cfs_dd2/checkSampleConsistencyVCF.py $MINOR \
#    -v ${SAMPLEPOST} -p FDK2D4 \
#    -s SM-7LV8O -s SM-7LV8P -s SM-7LV8Q  -s SM-7LV8R
#python $DISCO1/scripts/cfs_dd2/checkSampleConsistencyVCF.py $MINOR \
#    -v ${SAMPLEPOST} -p FDK2D4 \
#    -s SM-7LV8O -s SM-7LV8P -s SM-7LV8Q  -s SM-7LV8R
#if [[ $? != 0 ]]; then exit 1 ; fi

echo "CHECK SAMPLE CONSISTENCY DD2ALL PASS"
echo python $DISCO1/scripts/cfs_dd2/checkSampleConsistencyVCF.py $MINOR \
    -v ${SAMPLE/.vcf/.PASS.vcf} -p FDK2D4 \
    -s SM-7LV8O -s SM-7LV8P -s SM-7LV8Q  -s SM-7LV8R
python $DISCO1/scripts/cfs_dd2/checkSampleConsistencyVCF.py $MINOR \
    -v ${SAMPLE/.vcf/.PASS.vcf} -p FDK2D4 \
    -s SM-7LV8O -s SM-7LV8P -s SM-7LV8Q  -s SM-7LV8R
if [[ $? != 0 ]]; then exit 1 ; fi

echo "DD2FDK V REF (PFX) DD2 PASS ONLY" 
python $DISCO1/scripts/cfs_dd2/checkRealVFakeDD2.py \
    -v1 ${SAMPLE/.vcf/.PASS.FDK2D4.vcf} -v2 $COMPARISON $MINOR \
    -s SM-7LV8O -s SM-7LV8P -s SM-7LV8Q  \
    -p FDK \
    > ${SAMPLE/.vcf/.PASS.FDK2D4.FDKCONC.txt}

echo "DD22D4 V REF (PFX) DD2 PASS ONLY" 
python $DISCO1/scripts/cfs_dd2/checkRealVFakeDD2.py \
    -v1 ${SAMPLE/.vcf/.PASS.FDK2D4.FDKCONC.vcf} -v2 $COMPARISON $MINOR \
    -s SM-7LV8R \
    -p 2D4 \
    > ${SAMPLE/.vcf/.PASS.FDK2D4.FDKCONC.2D4CONC.txt}

#echo "REAL V FAKE DD2 ALL" 
#python $DISCO1/scripts/cfs_dd2/checkRealVFakeDD2.py \
#    -v1 ${SAMPLEPOST} -v2 $COMPARISON $MINOR \
#    -s SM-7LV8O -s SM-7LV8P -s SM-7LV8Q  -s SM-7LV8R \
#    > ${SAMPLE/.vcf/.DD2CONC.txt}


#cp ${SAMPLE/.vcf/.PASS.DD2CONC.vcf} ../${OUTFILE}.PASS.DD2CONC.vcf
cp ${SAMPLE/.vcf/.PASS.FDK2D4.FDKCONC.2D4CONC.vcf} ../${OUTFILE}.PASS.FDK2D4.FDKCONC.2D4CONC.vcf
cp  ${SAMPLE/.vcf/.failtab.comb.txt} ../${OUTFILE}.failtab.comb.txt
cp  ${SAMPLE/.vcf/.PASS.FDK2D4.FDKCONC.2D4CONC.txt} ../${OUTFILE}.DD2CONC.txt

# python $DISCO1/scripts/cfs/addSTRsToVCF.py -s ../STRs_mreps_Pf3D7_v3.R0.txt -v ../${OUTFILE}.PASS.DD2CONC.vcf

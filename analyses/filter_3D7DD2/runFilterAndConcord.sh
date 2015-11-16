#!/bin/bash

SAMPLE=$1
OUTFILE=$2
MINOR=$3

#minor = '-M' - passed directly to relevant scripts

SAMPLE_ABS=`readlink -f ${SAMPLE}`
SAMPLE=`basename $SAMPLE`


if [[ -z $OUTFILE ]] 
then
    OUTFILE=$SAMPLE
fi

mkdir tmp_${OUTFILE}_filter
cd tmp_${OUTFILE}_filter

ln -s $SAMPLE_ABS
ln -s $COMPARISON_ABS

REF="${WORK}/refs/Pf3D7_v3.fasta"
GATK=/humgen/gsa-hpprojects/GATK/bin/GenomeAnalysisTK-3.4-0-g7e26428

echo "FILTER HYPERVARIABLES"
python $DISCO1/scripts/cfs_dd2/filterHypervariableRecordsVCF.py -v ${SAMPLE}

SAMPLE=${SAMPLE/.gz/}

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

echo "CHECK SAMPLE CONSISTENCY 3D7Gold"
echo python $DISCO1/scripts/cfs_dd2/checkSampleConsistencyVCF.py \
    -v ${SAMPLE/.vcf/.HYPF.DD2LCHK.3D7LCHK.DD2FDK.vcf} -p 3D7gold  $MINOR \
    -s SM-7LV8S -s SM-7LV8T -s SM-7LV8U
python $DISCO1/scripts/cfs_dd2/checkSampleConsistencyVCF.py \
    -v ${SAMPLE/.vcf/.HYPF.DD2LCHK.3D7LCHK.DD2FDK.vcf} -p 3D7gold  $MINOR \
    -s SM-7LV8S -s SM-7LV8T -s SM-7LV8U
if [[ $? != 0 ]]; then exit 1 ; fi


SAMPLEPOST=${SAMPLE/.vcf/.HYPF.DD2LCHK.3D7LCHK.DD2FDK.3D7gold.vcf}

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


#echo "REAL V FAKE DD2 ALL" 
#python $DISCO1/scripts/cfs_dd2/checkRealVFakeDD2.py \
#    -v1 ${SAMPLEPOST} -v2 $COMPARISON $MINOR \
#    -s SM-7LV8O -s SM-7LV8P -s SM-7LV8Q  -s SM-7LV8R \
#    > ${SAMPLE/.vcf/.DD2CONC.txt}


#cp ${SAMPLE/.vcf/.PASS.DD2CONC.vcf} ../${OUTFILE}.PASS.DD2CONC.vcf

if [[ -f ${SAMPLE/.vcf/.PASS.vcf} ]]
then cp ${SAMPLE/.vcf/.PASS.vcf} ../${OUTFILE}.PASS.vcf
#elif  [[ -f ${SAMPLE/.vcf/.vcf} ]]
#then cp ${SAMPLE/.vcf/.PASS.vcf} ../${OUTFILE}.PASS.vcf
fi

cp  ${SAMPLE/.vcf/.failtab.comb.txt} ../${OUTFILE}.failtab.comb.txt
cp  ${SAMPLEPOST} ../${OUTFILE}.FILTERED.vcf

# python $DISCO1/scripts/cfs/addSTRsToVCF.py -s ../STRs_mreps_Pf3D7_v3.R0.txt -v ../${OUTFILE}.PASS.DD2CONC.vcf

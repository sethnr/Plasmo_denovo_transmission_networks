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

echo "CHECK LANE CONCORDANCE (all lanes)"
echo python $DISCO1/scripts/cfs_dd2/checkLaneConcordance.py \
    -v ${SAMPLE/.vcf/.HYPF.vcf} 
python $DISCO1/scripts/cfs_dd2/checkLaneConcordance.py \
    -v ${SAMPLE/.vcf/.HYPF.vcf} 
if [[ $? != 0 ]]; then exit 1 ; fi


SAMPLEPOST=${SAMPLE/.vcf/.HYPF.LCHK.vcf}

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

echo "CHECK SAMPLE CONSISTENCY (whole dataset)"
echo python $DISCO1/scripts/cfs_dd2/checkSampleConsistencyVCF.py $MINOR \
    -v ${SAMPLE/.vcf/.PASS.vcf} -p HFGR
python $DISCO1/scripts/cfs_dd2/checkSampleConsistencyVCF.py $MINOR \
    -v ${SAMPLE/.vcf/.PASS.vcf} -p HFGR 
if [[ $? != 0 ]]; then exit 1 ; fi


#cp ${SAMPLE/.vcf/.PASS.DD2CONC.vcf} ../${OUTFILE}.PASS.DD2CONC.vcf
cp ${SAMPLE/.vcf/.PASS.HFGR.vcf} ../${OUTFILE}.PASS.HFGR.vcf
cp  ${SAMPLE/.vcf/.failtab.comb.txt} ../${OUTFILE}.failtab.comb.txt
cp  ${SAMPLE/.vcf/.HYPF.LCHK.vcf} ../${OUTFILE}.FILTERS.vcf

# python $DISCO1/scripts/cfs/addSTRsToVCF.py -s ../STRs_mreps_Pf3D7_v3.R0.txt -v ../${OUTFILE}.PASS.DD2CONC.vcf
vcftools --vcf ${SAMPLE/.vcf/.PASS.HFGR.vcf} --max-missing 0.5 --recode --recode-INFO-all --out ${SAMPLE/.vcf/.HYPF.LCHK.miss0.5}

cp ${SAMPLE/.vcf/.HYPF.LCHK.miss0.5.recode.vcf} ../${OUTFILE}.PASS.HFGR.miss0.5.vcf


#!/bin/bash

REALVCF=$1
REFVCF=$2
PFXVCF=$3
OUTFILE=$4
MINOR=$5

#minor = '-M' - passed directly to relevant scripts

REALVCF_ABS=`readlink -f ${REALVCF}`
REALVCF=`basename $REALVCF`
REFVCF_ABS=`readlink -f ${REFVCF}`
REFVCF=`basename $REFVCF`
PFXVCF_ABS=`readlink -f ${PFXVCF}`
PFXVCF=`basename $PFXVCF`

if [[ -z $OUTFILE ]] 
then
    OUTFILE=$REALVCF
fi

mkdir tmp_${OUTFILE}_cfDD2
cd tmp_${OUTFILE}_cfDD2

ln -s $REALVCF_ABS
ln -s $REFVCF_ABS
ln -s $PFXVCF_ABS

REF="${WORK}/refs/Pf3D7_v3.fasta"
GATK=/humgen/gsa-hpprojects/GATK/bin/GenomeAnalysisTK-3.4-0-g7e26428

#####
# TECHNICAL FILTERS
#####

echo "FILTER NO-CALLS (100%)"
#vcftools --vcf $REALVCF -o ${REALVCF/.vcf/.NoCall0.vcf} --geno 1 --recode
python $DISCO1/scripts/cfs_dd2/filterNocallRecordsVCF.py -v $REALVCF
if [[ $? != 0 ]]; then exit 1 ; fi

echo "FILTER HYPERVARIABLES"
python $DISCO1/scripts/cfs_dd2/filterHypervariableRecordsVCF.py -v ${REALVCF/.vcf/.NoCall0.vcf}
if [[ $? != 0 ]]; then exit 1 ; fi

echo "CHECK LANE CONCORDANCE DD2"
echo python $DISCO1/scripts/cfs_dd2/checkLaneConcordance.py \
    -v ${REALVCF/.vcf/.NoCall0.HYPF.vcf} -p DD2 \
    -s SM-7LV8O -s SM-7LV8P -s SM-7LV8Q -s SM-7LV8R
python $DISCO1/scripts/cfs_dd2/checkLaneConcordance.py \
    -v ${REALVCF/.vcf/.NoCall0.HYPF.vcf} -p DD2 \
    -s SM-7LV8O -s SM-7LV8P -s SM-7LV8Q -s SM-7LV8R
if [[ $? != 0 ]]; then exit 1 ; fi

echo "CHECK LANE CONCORDANCE 3D7"
echo python $DISCO1/scripts/cfs_dd2/checkLaneConcordance.py \
    -v ${REALVCF/.vcf/.NoCall0.HYPF.DD2LCHK.vcf} -p 3D7 \
    -s SM-7LV8S -s SM-7LV8T -s SM-7LV8U -s SM-7LV8V
python $DISCO1/scripts/cfs_dd2/checkLaneConcordance.py \
    -v ${REALVCF/.vcf/.NoCall0.HYPF.DD2LCHK.vcf} -p 3D7 \
    -s SM-7LV8S -s SM-7LV8T -s SM-7LV8U -s SM-7LV8V
if [[ $? != 0 ]]; then exit 1 ; fi


#####
# BIO FILTERS
#####

echo "CHECK SAMPLE CONSISTENCY DD2FDK"
echo python $DISCO1/scripts/cfs_dd2/checkSampleConsistencyVCF.py \
    -v ${REALVCF/.vcf/.NoCall0.HYPF.DD2LCHK.3D7LCHK.vcf} -p DD2FDK  $MINOR \
    -s SM-7LV8O -s SM-7LV8P -s SM-7LV8Q
python $DISCO1/scripts/cfs_dd2/checkSampleConsistencyVCF.py \
    -v ${REALVCF/.vcf/.NoCall0.HYPF.DD2LCHK.3D7LCHK.vcf} -p DD2FDK  $MINOR \
    -s SM-7LV8O -s SM-7LV8P -s SM-7LV8Q
if [[ $? != 0 ]]; then exit 1 ; fi

echo "CHECK SAMPLE CONSISTENCY 3D7Gold"
echo python $DISCO1/scripts/cfs_dd2/checkSampleConsistencyVCF.py \
    -v ${REALVCF/.vcf/.NoCall0.HYPF.DD2LCHK.3D7LCHK.DD2FDK.vcf} -p 3D7GOLD  $MINOR \
    -s SM-7LV8S -s SM-7LV8PT -s SM-7LV8U
python $DISCO1/scripts/cfs_dd2/checkSampleConsistencyVCF.py \
    -v ${REALVCF/.vcf/.NoCall0.HYPF.DD2LCHK.3D7LCHK.DD2FDK.vcf} -p 3D7GOLD  $MINOR \
    -s SM-7LV8S -s SM-7LV8T -s SM-7LV8U
if [[ $? != 0 ]]; then exit 1 ; fi


REALVCFPOST=${REALVCF/.vcf/.NoCall0.HYPF.DD2LCHK.3D7LCHK.DD2FDK.3D7GOLD.vcf}
echo "FILTERING:"
perl $DISCO1/scripts/vcf/countFilters.pl ${REALVCFPOST} >  ${REALVCF/.vcf/.failtab.comb.txt}
java -jar $GATK/GenomeAnalysisTK.jar \
   -T SelectVariants \
   -R $REF \
   -V ${REALVCFPOST} \
   -o ${REALVCF/.vcf/.PASS.vcf} \
   -env \
   -ef \
   -trimAlternates

REALVCF=${REALVCF/.vcf/.PASS.vcf}
########
# WITHIN FILTERED SET, LOOK AT DISCORDANTS
########
# NB: have already removed internal discords in goldberg & Fidock reps
# hence any new discords are 2D4vFDK or Goldberg V Wirth

echo "CHECK SAMPLE CONSISTENCY DD2: 2D4 v FDK"
echo python $DISCO1/scripts/cfs_dd2/checkSampleConsistencyVCF.py $MINOR \
    -v $REALVCF -p 2D4FDK \
    -s SM-7LV8O -s SM-7LV8P -s SM-7LV8Q  -s SM-7LV8R
python $DISCO1/scripts/cfs_dd2/checkSampleConsistencyVCF.py $MINOR \
    -v $REALVCF -p 2D4FDK \
    -s SM-7LV8O -s SM-7LV8P -s SM-7LV8Q  -s SM-7LV8R
if [[ $? != 0 ]]; then exit 1 ; fi

echo "CHECK SAMPLE CONSISTENCY 3D7: Goldberg v Wirth"
echo python $DISCO1/scripts/cfs_dd2/checkSampleConsistencyVCF.py $MINOR \
    -v ${REALVCF/.vcf/.2D4FDK.vcf} -p GoldWirth \
    -s SM-7LV8S -s SM-7LV8T -s SM-7LV8U  -s SM-7LV8V
python $DISCO1/scripts/cfs_dd2/checkSampleConsistencyVCF.py $MINOR \
    -v ${REALVCF/.vcf/.2D4FDK.vcf} -p GoldWirth \
    -s SM-7LV8S -s SM-7LV8T -s SM-7LV8U  -s SM-7LV8V
if [[ $? != 0 ]]; then exit 1 ; fi

#VS REF DISCOVAR
echo "REAL V REF DD2: Fidock" 
python $DISCO1/scripts/cfs_dd2/checkRealVFakeDD2.py \
    -v1 ${REALVCF/.vcf/.2D4FDK.GoldWirth.vcf} -v2 $REFVCF $MINOR \
    -s SM-7LV8O -s SM-7LV8P -s SM-7LV8Q  \
    -p Fdk \
    > /dev/null
if [[ $? != 0 ]]; then exit 1 ; fi

echo "REAL V REF DD2: 2D4" 
python $DISCO1/scripts/cfs_dd2/checkRealVFakeDD2.py \
    -v1 ${REALVCF/.vcf/.2D4FDK.GoldWirth.FdkCONC.vcf} -v2 $REFVCF $MINOR \
    -s SM-7LV8R \
    -p 2D4 \
    > /dev/null
if [[ $? != 0 ]]; then exit 1 ; fi

#VS PFX
echo "REAL V PFX DD2: Fidock" 
python $DISCO1/scripts/cfs_dd2/checkRealVFakeDD2.py \
    -v1 ${REALVCF/.vcf/.2D4FDK.GoldWirth.FdkCONC.2D4CONC.vcf} -v2 $PFXVCF $MINOR \
    -s SM-7LV8O -s SM-7LV8P -s SM-7LV8Q  \
    -p FdkX \
    > /dev/null
if [[ $? != 0 ]]; then exit 1 ; fi

echo "REAL V PFX DD2: 2D4" 
python $DISCO1/scripts/cfs_dd2/checkRealVFakeDD2.py \
    -v1 ${REALVCF/.vcf/.2D4FDK.GoldWirth.FdkCONC.2D4CONC.FdkXCONC.vcf} -v2 $PFXVCF $MINOR \
    -s SM-7LV8R \
    -p 2D4X \
    > ${REALVCF/.vcf/.RefPFXConc.txt}
if [[ $? != 0 ]]; then exit 1 ; fi

REALVCFPOST=${REALVCF/.vcf/.2D4FDK.GoldWirth.FdkCONC.2D4CONC.FdkXCONC.2D4XCONC.vcf} 
perl $DISCO1/scripts/vcf/countFilters.pl ${REALVCFPOST} >  ${REALVCF/.vcf/.discord.comb.txt}


cp ${REALVCF/.vcf/.2D4FDK.GoldWirth.FdkCONC.2D4CONC.FdkXCONC.2D4XCONC.vcf} ../${OUTFILE}.PASS.RefPFXConc.vcf
cp  ${REALVCF/.PASS.vcf/.failtab.comb.txt} ../${OUTFILE}.failtab.comb.txt
cp  ${REALVCF/.vcf/.discord.comb.txt} ../${OUTFILE}.discord.comb.txt
cp  ${REALVCF/.vcf/.RefPFXConc.txt} ../${OUTFILE}.RefPFXConc.txt



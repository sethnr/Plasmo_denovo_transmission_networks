#!/bin/bash

FAKEDD2=$1
REALDD2=$2
OUTFILE=$3
MINOR=$4

REALDD2_ABS=`readlink -f ${REALDD2}`
FAKEDD2_ABS=`readlink -f ${FAKEDD2}`
REALDD2=`basename $REALDD2`
FAKEDD2=`basename $FAKEDD2`

if [[ -z $OUTFILE ]] 
then
    OUTFILE=$REALDD2
fi

mkdir tmp_${OUTFILE}_cfDD2
cd tmp_${OUTFILE}_cfDD2

ln -s $REALDD2_ABS
ln -s $FAKEDD2_ABS

REF="${WORK}/refs/Pf3D7_v3.fasta"
GATK=/humgen/gsa-hpprojects/GATK/bin/GenomeAnalysisTK-3.4-0-g7e26428

echo "FILTER HYPERVARIABLES"
python $DISCO1/scripts/cfs_dd2/filterHypervariableRecordsVCF.py -v ${REALDD2}

echo "CHECK LANE CONCORDANCE DD2"
python $DISCO1/scripts/cfs_dd2/checkLaneConcordance.py \
    -v ${REALDD2/.vcf/.HYPF.vcf} -p DD2 \
    -s SM-7LV8O -s SM-7LV8P -s SM-7LV8Q -s SM-7LV8R
if [[ $? != 0 ]]; then exit 1 ; fi

echo "CHECK LANE CONCORDANCE 3D7"
python $DISCO1/scripts/cfs_dd2/checkLaneConcordance.py \
    -v ${REALDD2/.vcf/.HYPF.DD2LCHK.vcf} -p 3D7 \
    -s SM-7LV8S -s SM-7LV8T -s SM-7LV8U -s SM-7LV8V
if [[ $? != 0 ]]; then exit 1 ; fi

echo "CHECK SAMPLE CONSISTENCY DD2FDK"
python $DISCO1/scripts/cfs_dd2/checkSampleConsistencyVCF.py \
    -v ${REALDD2/.vcf/.HYPF.DD2LCHK.3D7LCHK.vcf} -p DD2FDK  $MINOR \
    -s SM-7LV8O -s SM-7LV8P -s SM-7LV8Q
if [[ $? != 0 ]]; then exit 1 ; fi

echo "CHECK SAMPLE CONSISTENCY DD2ALL"
python $DISCO1/scripts/cfs_dd2/checkSampleConsistencyVCF.py $MINOR \
    -v ${REALDD2/.vcf/.HYPF.DD2LCHK.3D7LCHK.DD2FDK.vcf} -p DD2ALL \
    -s SM-7LV8O -s SM-7LV8P -s SM-7LV8Q  -s SM-7LV8R
if [[ $? != 0 ]]; then exit 1 ; fi

#python $DISCO1/scripts/cfs_dd2/filterHypervariableRecordsVCF.py -v ${REALDD2/.vcf/.DD2LCHK.3D7LCHK.DD2FDK.DD2ALL.vcf}


REALDD2POST=${REALDD2/.vcf/.HYPF.DD2LCHK.3D7LCHK.DD2FDK.DD2ALL.vcf}

#count filters
echo "FILTERING:"
#perl -lane 'BEGIN{my %filt;} next if (m/^\#/); $filt{$F[6]}++; END{foreach $f (keys %filt) {print $f."\t".$filt{$f} }}' ${REALDD2POST}   | sort > ${REALDD2}.failtab.comb.txt

perl $DISCO1/scripts/vcf/countFilters.pl ${REALDD2POST} >  ${REALDD2/.vcf/.failtab.comb.txt}

#remove filtered and non-varying loci
#remove non-varying alleles (hypervariables?)
java -jar $GATK/GenomeAnalysisTK.jar \
   -T SelectVariants \
   -R $REF \
   -V ${REALDD2POST} \
   -o ${REALDD2/.vcf/.PASS.vcf} \
   -env \
   -ef \
   -trimAlternates

echo "REAL V FAKE DD2 PASS ONLY" 
python $DISCO1/scripts/cfs_dd2/checkRealVFakeDD2.py \
    -v1 ${REALDD2/.vcf/.PASS.vcf} -v2 $FAKEDD2 $MINOR \
    -s SM-7LV8O -s SM-7LV8P -s SM-7LV8Q  -s SM-7LV8R \
    > ${REALDD2/.vcf/.PASS.DD2CONC.txt}

echo "REAL V FAKE DD2 ALL" 
python $DISCO1/scripts/cfs_dd2/checkRealVFakeDD2.py \
    -v1 ${REALDD2POST} -v2 $FAKEDD2 $MINOR \
    -s SM-7LV8O -s SM-7LV8P -s SM-7LV8Q  -s SM-7LV8R \
    > ${REALDD2/.vcf/.DD2CONC.txt}


cp ${REALDD2/.vcf/.PASS.DD2CONC.vcf} ../${OUTFILE}.PASS.DD2CONC.vcf
cp  ${REALDD2/.vcf/.failtab.comb.txt} ../${OUTFILE}.failtab.comb.txt
cp  ${REALDD2/.vcf/.DD2CONC.txt} ../${OUTFILE}.DD2CONC.txt

#!/bin/bash

VCF=$1
STR=$2
if [[ -z "$STR" ]]; then 
  STR="$DISCO1/analyses/STR_calling/STRs_mreps_Pf3D7_v3.R0.txt"
fi

CONSDB=$3
OUTFILE=$4
if [[ -z "$CONSDB" ]]; then 
  CONSDB="$DISCO1/analyses/SNPconsequences/SNPeff/Pf3D7_v3_a25"
fi

VCF_ABS=`readlink -f ${VCF}`
STR_ABS=`readlink -f ${STR}`
CONS_ABS=`readlink -f ${CONSDB}`
VCF=`basename $VCF`
STR=`basename $STR`
CONSDB=`basename $CONSDB`

if [[ -z $OUTFILE ]] 
then
    OUTFILE=${VCF/.vcf/.VARCLASS}
fi

mkdir tmp_${OUTFILE}_varClass
cd tmp_${OUTFILE}_varClass

ln -s $VCF_ABS
ln -s $STR_ABS
ln -s $CONS_ABS

REF="${WORK}/refs/Pf3D7_v3.fasta"
GATK=/humgen/gsa-hpprojects/GATK/bin/GenomeAnalysisTK-3.4-0-g7e26428

echo "ADDING STR INFO TO VCF"
python $DISCO1/scripts/cfs/addSTRsToVCF.py -s $STR -v $VCF

echo "ADDING GENE/CONS INFO TO VCF"
snpEff eff $CONSDB ${VCF/.vcf/.STRs.vcf} \
       -no-downstream -no-upstream -no-intron -no-intergenic \
       > ${VCF/.vcf/.STRs.ANN.vcf}
cp ${VCF/.vcf/.STRs.ANN.vcf} ../${OUTFILE}.STRs.ANN.vcf


# echo "BLOCK SUMMARISE VCF"
# python $DISCO1/scripts/vcf/blockSummariseVCF.py \
#    -i ANN:Count -i STR:Count -i DD2ALLMAF:Mean \
#    -v ${VCF/.vcf/.STRs.ANN.vcf}

echo "SUMMARISE VAR CLASSES IN VCF"
python $DISCO1/scripts/vcf/summariseVarClassesVCF.py \
       -v ${VCF/.vcf/.STRs.ANN.vcf} > ${VCF/.vcf/.varClass.txt}

cp  ${VCF/.vcf/.varClass.txt} ../${OUTFILE}.varClass.txt
#cp  ${REALDD2/.vcf/.DD2CONC.txt} ../${OUTFILE}.DD2CONC.txt


#!/bin/bash

OUT="BEASTer"

while getopts "d:R:l:o:V:D" opt; do
  case $opt in
    l) LOCUS=$OPTARG;;
    R) REF=${OPTARG};;
    o) OUT=$OPTARG;;
    V) VCF=$OPTARG;;
    D) DELETE=1;;
   \?) echo "Invalid option: -$OPTARG" >&2 ;;
  esac
done
shift $(expr $OPTIND - 1 )

#if array job, make TMP folder for each
if [ -n "$LSB_JOBINDEX" ]; then
    JOBINDEX=${LSB_JOBINDEX}
    OUT=${OUT}.${JOBINDEX}
elif [ -n "$SGE_TASK_ID" ]; then
    JOBINDEX=${SGE_TASK_ID}
    OUT=${OUT}.${JOBINDEX}
else
    OUT=${OUT}
fi


echo $REF

#echo samtools faidx $REF $LOCUS \> ${LOCUS}.fasta
#samtools faidx $REF $LOCUS > ${LOCUS}.fasta
#tabix -h $VCF $LOCUS > ${LOCUS}.vcf
#bgzip ${LOCUS}.vcf
#tabix -pvcf ${LOCUS}.vcf.gz
#splitVcfBySample.sh ${LOCUS}.vcf.gz $REF
#samtools faidx $REF $LOCUS > locus.fasta
#tabix -h $VCF $LOCUS > locus.vcf
#bgzip locus.vcf
#tabix -pvcf locus.vcf.gz
#splitVcfBySample.sh locus.vcf.gz $REF

if [[ "$#" > 0 ]]
then
    echo "passed samples"
    SAMPLES=$@
else

    if [[ $VCF == *".gz" ]]
    then 
	SAMPLES=`zcat $VCF | head -n 1000 | grep \#CHR | cut -f 10- `
    else
	SAMPLES=`head -n 1000 ${VCF} | grep \#CHR | cut -f 10- `
    fi
fi

VCF=`readlink -f $VCF`

mkdir tmp_${OUT}
cd tmp_${OUT}


VARCOUNT=`tabix $VCF $LOCUS | wc -l `
tabix -h $VCF $LOCUS > locus.vcf 
bgzip locus.vcf 
tabix -pvcf locus.vcf.gz 

rm locus.fasta
for SAMP in $SAMPLES ; 
do  

  echo samtools faidx $REF $LOCUS \| vcf-consensus -s $SAMP locus.vcf.gz \> ${SAMP}.fasta
  samtools faidx $REF $LOCUS | vcf-consensus -s $SAMP locus.vcf.gz 1> ${SAMP}.fasta  & wait
  perl -i -pe "s/\>.*\n/\>${SAMP}\n/" ${SAMP}.fasta & wait
  cat ${SAMP}.fasta >> locus.fasta
done

wc -l *fasta


echo "clustalw2 -infile=locus.fasta -outfile=locus.nexus -output=nexus"
clustalw2 -infile=locus.fasta -outfile=locus.nexus -output=nexus & wait

TEMPLATE=/seq/plasmodium/sredmond/pfdisco/analyses/substitution_rates/substitution_rate.template
#echo "./BEASTGen/bin/beastgen -date_regex '\d\d$' BEASTGen/templates/substitution_rate.template locus.nexus locus.xml"
#~/bin/beastgen 
#beastgen -date_regex '\d\d$' BEASTGen/templates/substitution_rate.template locus.nexus locus.xml
cp $TEMPLATE ./
echo beastgen -date_regex '\d\d$' substitution_rate.template  locus.nexus locus.xml
beastgen -date_regex '\d\d$' substitution_rate.template  locus.nexus locus.xml & wait

echo "java -jar BEASTv1.8.3/lib/beast.jar locus.xml \> locus.beast.txt"
#java -jar ~/software/BEASTv1.8.3/lib/beast.jar locus.xml > locus.beast.txt
beast -threads 3 locus.xml > locus.beast.txt
#beast locus.xml & wait


MEANRATE=`perl -lane 'BEGIN{$C=0; $T=0} next if $F[0] =~ m/\#/gi || $F[0]<1e5; $C+=1; $T+=$F[11]; END{print ($T/$C)}' locus.log`
#echo perl -lane \'$MEANRATE\' locus.log
#perl -lane 'BEGIN{$C=0; $T=0} next if $F[0] =~ m/\#/gi || $F[0]<1e5; $C+=1; $T+=$F[11]; END{print join("\t",${LOCUS},($T/$C))}' locus.log
#perl -lane \' $MEANRATE \' locus.log
echo "$OUT    $LOCUS    $VARCOUNT    $MEANRATE" > ../${OUT}_beastER
echo "$OUT    $LOCUS    $VARCOUNT    $MEANRATE" > ${OUT}_beastER

cd ../

if [[ -n "$DELETE" ]]
then
  echo "DELETING tmp_${OUT}"    
  rm -r tmp_${OUT}    
fi

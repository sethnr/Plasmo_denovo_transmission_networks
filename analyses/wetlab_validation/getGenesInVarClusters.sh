#!/bin/bash

VCF=$1
GFF=$2

MERGE=10000
FLANK=1000

BED=${VCF/.vcf/.REGIONS.bed}

bedtools merge -i $VCF -d $MERGE -n > $BED
#echo "vars"
#head $BED
perl -i -lane "\$F[1] = \$F[1]-${FLANK}; \$F[2] = \$F[2]+${FLANK} ; print join(\"\t\",@F);" $BED
#echo "regions"
#head $BED

bedtools intersect -a $GFF -b $BED  2>&1> ${VCF/.vcf/.genes}

perl -i -lane 'BEGIN{use URI::Escape;}; if ($F[1]=~ m/^\d+$/) {print $_} else {$F[8] =~ s/.*Name=([^;\s]+).*description=([^;\s]+).*/$1\t$2/gi; $F[8] =~ s/\+/\ /gi; print join("\t",$F[0],$F[3],$F[4],uri_unescape($F[8]))}'  ${VCF/.vcf/.genes}

cat $BED ${VCF/.vcf/.genes} | sort -k1,1 -k2,2g
#sort -k2,2g ${VCF/.vcf/.intersect.s}


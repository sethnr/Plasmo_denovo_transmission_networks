#!/bin/bash

VCF=$1

python checkTemporalConsistencyVCF.py -v $VCF \
  -s DD2-2D4 -s HFGRIII.1-30x -s HFGRIII.1-60x \
  -s HFGRIII.2-60x -s HFGRIII.3-200x > ${VCF/.vcf.gz/.HFIII.temporal.txt}

python checkTemporalConsistencyVCF.py -v $VCF \
  -s DD2-2D4  -s HFGRII.2-10x -s HFGRII.1-30x  \
  -s HFGRII.1-200x -s HFGRII.2-200x \
  -s HFGRII.2-300x  > ${VCF/.vcf.gz/.HFII.temporal.txt}

perl -i -lane 'print $_ if $F[5] eq "0" && $F[12] > 0.5' ${VCF/.vcf.gz/.HFII.temporal.txt}
perl -i -lane 'print $_ if $F[5] eq "0" && $F[12] > 0.5' ${VCF/.vcf.gz/.HFIII.temporal.txt}

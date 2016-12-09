#!/bin/bash

perl -lane 'print join(" ","ln -s",$F[4],$F[0]."_".$F[1].".bam"); $F[4] =~ s/bam/bai/; print join(" ","ln -s",$F[4],$F[0]."_".$F[1].".bam.bai");' /seq/plasmodium/sredmond/pfdisco/samples/3D7DD2.samples.txt > getAllSymlinks.sh



for SAMPLE in `cut -f 1 /seq/plasmodium/sredmond/pfdisco/samples/HFGR.samples.txt | uniq`;
do
    bsub13 cnator_${SAMPLE}  $WORK/idi_broad_scripts/wrappers/runCNVNATOR.sh -o $SAMPLE data/${SAMPLE}*bam
done

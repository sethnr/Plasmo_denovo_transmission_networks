#!/bin/bash


FPS=""
FCS=""
PCS=""
for SAMPLE in `cut -f 1 /seq/plasmodium/sredmond/pfdisco/samples/HFGR.samples.txt | uniq`
do 
    python mergeCNVcalls.py --mean \
       -f freec/${SAMPLE}_CNVs \
       -c cnvnator/${SAMPLE}_CNVs \
       > comb_${SAMPLE}_FC.CNVs
    FCS="${FCS} -f comb_${SAMPLE}_FC.CNVs"
done
python mergeCNVSamples.py $FCS > comb_FCs.CNVs.tab


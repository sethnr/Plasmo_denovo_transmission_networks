#!/bin/bash 

#GET RELEVANT FRAGMENT
samtools faidx $WORK/refs/PfIT_v3.fasta PfIT_07_v3:450000-750000 > PfIT_07_v3:450000-750000.fasta
samtools faidx $WORK/refs/Pf3D7_v3.fasta Pf3D7_07_v3:450000-750000 > Pf3D7_07_v3:450000-750000.fasta
samtools faidx $WORK/refs/PfDD2_v1.fasta  PfDd2_07_TT:450000-750000 > PfDD2_07_v1:450000-750000.fasta

bgzip PfIT_07_v3:450000-750000.fasta
bgzip Pf3D7_07_v3:450000-750000.fasta
bgzip PfDD2_07_v1:450000-750000.fasta

#WGSIM
$DISCO1/scripts/makeSimSeq.sh -q Pf3D7_07_v3:450000-750000.fasta.gz
$DISCO1/scripts/makeSimSeq.sh -q PfDD2_07_v1:450000-750000.fasta.gz
$DISCO1/scripts/makeSimSeq.sh -q PfIT_07_v3:450000-750000.fasta.gz

#ALIGN ALL
bsub13 alsim.DD23D7 \
$DISCO1/scripts/alignSimSeq.sh -q PfDD2_07_v1:450000-750000.fasta.gz -r $WORK/refs/Pf3D7_v3.fasta -o PfDD2_07_v1:450000-750000_v_Pf3D7
bsub13 alsim.DD2DD2 \
$DISCO1/scripts/alignSimSeq.sh -q PfDD2_07_v1:450000-750000.fasta.gz -r $WORK/refs/PfDD2_v1.fasta -o PfDD2_07_v1:450000-750000_v_PfDD2
bsub13 alsim.DD2IT  \
$DISCO1/scripts/alignSimSeq.sh -q PfDD2_07_v1:450000-750000.fasta.gz -r $WORK/refs/PfIT_v3.fasta  -o PfDD2_07_v1:450000-750000_v_PfIT

bsub13 alsim.3D73D7   \
$DISCO1/scripts/alignSimSeq.sh -q Pf3D7_07_v3:450000-750000.fasta.gz -r $WORK/refs/Pf3D7_v3.fasta -o Pf3D7_07_v3:450000-750000_v_Pf3D7
bsub13 alsim.3D7DD2  \
$DISCO1/scripts/alignSimSeq.sh -q Pf3D7_07_v3:450000-750000.fasta.gz -r $WORK/refs/PfDD2_v1.fasta -o Pf3D7_07_v3:450000-750000_v_PfDD2
bsub13 alsim.3D7IT  \
$DISCO1/scripts/alignSimSeq.sh -q Pf3D7_07_v3:450000-750000.fasta.gz -r $WORK/refs/PfIT_v3.fasta  -o Pf3D7_07_v3:450000-750000_v_PfIT

bsub13 alsim.IT3D7 \
$DISCO1/scripts/alignSimSeq.sh -q PfIT_07_v3:450000-750000.fasta.gz -r $WORK/refs/Pf3D7_v3.fasta -o PfIT_07_v3:450000-750000_v_Pf3D7
bsub13 alsim.ITDD2 \
$DISCO1/scripts/alignSimSeq.sh -q PfIT_07_v3:450000-750000.fasta.gz -r $WORK/refs/PfDD2_v1.fasta -o PfIT_07_v3:450000-750000_v_PfDD2
bsub13 alsim.ITIT \
$DISCO1/scripts/alignSimSeq.sh -q PfIT_07_v3:450000-750000.fasta.gz -r $WORK/refs/PfIT_v3.fasta  -o PfIT_07_v3:450000-750000_v_PfIT


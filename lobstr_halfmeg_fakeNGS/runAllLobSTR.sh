bsub13 lobSTR.r3D7.q3D7 $DISCO1/scripts/runLobSTR.sh \
   -b $DISCO1/analyses/fakeNGS/Pf3D7_07_v3:450000-750000_v_Pf3D7.bam \
   -o lobSTR.r3D7.q3D7

bsub13 lobSTR.rDD2.qDD2 $DISCO1/scripts/runLobSTR.sh \
   -b $DISCO1/analyses/fakeNGS/PfDD2_07_v1:450000-750000_v_PfDD2.bam \
   -o lobSTR.rDD2.qDD2 \
   -r $WORK/lobSTRrefs/PfDD2_lobstr


bsub13 lobSTR.r3D7.qDD2 $DISCO1/scripts/runLobSTR.sh \
   -b $DISCO1/analyses/fakeNGS/Pf3D7_07_v3:450000-750000_v_PfDD2.bam \
   -o lobSTR.r3D7.qDD2 -s PfDD2
bsub13 lobSTR.rDD2.q3D7 $DISCO1/scripts/runLobSTR.sh \
   -b $DISCO1/analyses/fakeNGS/PfDD2_07_v1:450000-750000_v_Pf3D7.bam \
   -o lobSTR.rDD2.q3D7 -s Pf3D7 -l Pf3D7\
   -r $WORK/lobSTRrefs/PfDD2_lobstr

addReadGroups -s Pf3D7 -s Pf3D7 -l Pf3D7 -b Pf3D7_07_v3:450000-750000_v_Pf3D7.bam -o Pf3D7_07_v3:450000-750000_v_Pf3D7.rg.bam -i H2MGC.1 -I0
addReadGroups -s Pf3D7 -s Pf3D7 -l Pf3D7 -b Pf3D7_07_v3:450000-750000_v_PfDD2.bam -o Pf3D7_07_v3:450000-750000_v_PfDD2.rg.bam -i H2MGC.1 -I0
addReadGroups -s Pf3D7 -s Pf3D7 -l Pf3D7 -b Pf3D7_07_v3:450000-750000_v_PfIT.bam -o Pf3D7_07_v3:450000-750000_v_PfIT.rg.bam -i H2MGC.1 -I0

addReadGroups -s PfDD2 -s PfDD2 -l PfDD2 -b PfDD2_07_v1:450000-750000_v_Pf3D7.bam -o PfDD2_07_v1:450000-750000_v_Pf3D7.rg.bam -i H2MGC.1 -I0
addReadGroups -s PfDD2 -s PfDD2 -l PfDD2 -b PfDD2_07_v1:450000-750000_v_PfDD2.bam -o PfDD2_07_v1:450000-750000_v_PfDD2.rg.bam -i H2MGC.1 -I0
addReadGroups -s PfDD2 -s PfDD2 -l PfDD2 -b PfDD2_07_v1:450000-750000_v_PfIT.bam -o PfDD2_07_v1:450000-750000_v_PfIT.rg.bam -i H2MGC.1 -I0

addReadGroups -s PfIT -s PfIT -l PfIT -b PfIT_07_v3:450000-750000_v_Pf3D7.bam -o PfIT_07_v3:450000-750000_v_Pf3D7.rg.bam -i H2MGC.1 -I0
addReadGroups -s PfIT -s PfIT -l PfIT -b PfIT_07_v3:450000-750000_v_PfDD2.bam -o PfIT_07_v3:450000-750000_v_PfDD2.rg.bam -i H2MGC.1 -I0
addReadGroups -s PfIT -s PfIT -l PfIT -b PfIT_07_v3:450000-750000_v_PfIT.bam -o PfIT_07_v3:450000-750000_v_PfIT.rg.bam -i H2MGC.1 -I0

for BAM in *rg.bam; do samtools sort -o ${BAM/.bam/.s.bam} -T temp  $BAM; done
for BAM in *s.bam; do mv $BAM ${BAM/.s.bam/.bam}; done
for BAM in *bam; do samtools index $BAM; done

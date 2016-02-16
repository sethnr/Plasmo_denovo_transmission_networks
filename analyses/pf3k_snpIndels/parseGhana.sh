for CHR in 01 02 03 04 05 06 07 08 09 10 11 12 13 14;
do 
  vcftools --gzvcf SNP_INDEL_Pf3D7_${CHR}_v3.combined.filtered.vcf.gz \
	--thin 1000 --mac 1 \
	--keep $DISCO1/analyses/pf3k_snpIndels/Pf3k_ghana_ids.txt \
	--out SNP_INDEL_Pf3D7_${CHR}_v3_ghana_1K \
	--recode
  mv SNP_INDEL_Pf3D7_${CHR}_v3_ghana_1K.recode.vcf SNP_INDEL_Pf3D7_${CHR}_v3_ghana_1K.vcf 
  bgzip SNP_INDEL_Pf3D7_${CHR}_v3_ghana_1K.vcf
  tabix -pvcf SNP_INDEL_Pf3D7_${CHR}_v3_ghana_1K.vcf.gz
done

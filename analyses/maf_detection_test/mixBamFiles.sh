#!/bin/bash


ln -s /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/1/SM-7LV8R/H2MGCBCXX.1.aligned.duplicates_marked.bam Dd22D4_1.bam
ln -s /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/2/SM-7LV8R/H2MGCBCXX.2.aligned.duplicates_marked.bam Dd22D4_2.bam
ln -s /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/1/SM-7LV8P/H2MGCBCXX.1.aligned.duplicates_marked.bam Dd2FDK2_1.bam
ln -s /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/2/SM-7LV8P/H2MGCBCXX.2.aligned.duplicates_marked.bam Dd2FDK2_2.bam

ln -s /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/1/SM-7LV8R/H2MGCBCXX.1.aligned.duplicates_marked.bai Dd22D4_1.bai
ln -s /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/2/SM-7LV8R/H2MGCBCXX.2.aligned.duplicates_marked.bai Dd22D4_2.bai
ln -s /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/1/SM-7LV8P/H2MGCBCXX.1.aligned.duplicates_marked.bai Dd2FDK2_1.bai
ln -s /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/2/SM-7LV8P/H2MGCBCXX.2.aligned.duplicates_marked.bai Dd2FDK2_2.bai


mkdir data

for BAM in *.bam
do
    for FRACT in 0.02 0.05 0.1 0.2 0.5 0.8 0.9 0.95 0.98 1
    do
	echo $FRACT, $BAM
	samtools view -s $FRACT -b $BAM Pf3D7_04_v3 > data/${BAM/.bam/.c4.$FRACT.bam}
done
done


mkdir comb
echo "MERGING 0.02"
samtools merge -r comb/Dd2comb_2D4_0.02.bam  data/Dd22D4_1.c4.0.02.bam data/Dd2FDK2_1.c4.0.98.bam \
    data/Dd22D4_2.c4.0.02.bam data/Dd2FDK2_2.c4.0.98.bam

echo "MERGING 0.05"
samtools merge -r comb/Dd2comb_2D4_0.05.bam  data/Dd22D4_1.c4.0.05.bam data/Dd2FDK2_1.c4.0.95.bam \
    data/Dd22D4_2.c4.0.05.bam data/Dd2FDK2_2.c4.0.95.bam
echo "MERGING 0.1"
samtools merge -r comb/Dd2comb_2D4_0.1.bam  data/Dd22D4_1.c4.0.2.bam data/Dd2FDK2_1.c4.0.9.bam \
    data/Dd22D4_2.c4.0.2.bam data/Dd2FDK2_2.c4.0.9.bam
echo "MERGING 0.2"
samtools merge -r comb/Dd2comb_2D4_0.2.bam  data/Dd22D4_1.c4.0.2.bam data/Dd2FDK2_1.c4.0.8.bam \
    data/Dd22D4_2.c4.0.2.bam data/Dd2FDK2_2.c4.0.8.bam
echo "MERGING 0.5"
samtools merge -r comb/Dd2comb_2D4_0.5.bam  data/Dd22D4_1.c4.0.5.bam data/Dd2FDK2_1.c4.0.5.bam \
    data/Dd22D4_2.c4.0.5.bam data/Dd2FDK2_2.c4.0.5.bam
echo "MERGING 1.0"
samtools merge -r comb/Dd2comb_2D4_1.0.bam  data/Dd22D4_1.c4.1.bam data/Dd22D4_2.c4.1.bam



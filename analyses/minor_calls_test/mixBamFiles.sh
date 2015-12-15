#!/bin/bash


ln -s /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/1/SM-7LV8R/H2MGCBCXX.1.aligned.duplicates_marked.bam Dd22D4_1.bam
ln -s /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/2/SM-7LV8R/H2MGCBCXX.2.aligned.duplicates_marked.bam Dd22D4_2.bam
ln -s /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/1/SM-7LV8P/H2MGCBCXX.1.aligned.duplicates_marked.bam Dd2FDK2_1.bam
ln -s /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/2/SM-7LV8P/H2MGCBCXX.2.aligned.duplicates_marked.bam Dd2FDK2_2.bam

ln -s /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/1/SM-7LV8R/H2MGCBCXX.1.aligned.duplicates_marked.bai Dd22D4_1.bai
ln -s /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/2/SM-7LV8R/H2MGCBCXX.2.aligned.duplicates_marked.bai Dd22D4_2.bai
ln -s /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/1/SM-7LV8P/H2MGCBCXX.1.aligned.duplicates_marked.bai Dd2FDK2_1.bai
ln -s /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/2/SM-7LV8P/H2MGCBCXX.2.aligned.duplicates_marked.bai Dd2FDK2_2.bai


for BAM in *.bam;
do
    for FRACT in 0.2 0.4 0.5 0.6 0.8;
    do
#	samtools view -s $FRACT -b $BAM Pf3D7_04_v3 > data/${BAM/.bam/.c4.$FRACT.bam}
    done
done


samtools cat -o data/Dd2comb_2D40.2_L1.bam  data/Dd22D4_1.c4.0.2.bam data/Dd2FDK2_1.c4.0.8.bam
samtools cat -o data/Dd2comb_2D40.4_L1.bam  data/Dd22D4_1.c4.0.4.bam data/Dd2FDK2_1.c4.0.6.bam
samtools cat -o data/Dd2comb_2D40.5_L1.bam  data/Dd22D4_1.c4.0.5.bam data/Dd2FDK2_1.c4.0.5.bam
samtools cat -o data/Dd2comb_2D40.6_L1.bam  data/Dd22D4_1.c4.0.6.bam data/Dd2FDK2_1.c4.0.4.bam
samtools cat -o data/Dd2comb_2D40.8_L1.bam  data/Dd22D4_1.c4.0.8.bam data/Dd2FDK2_1.c4.0.2.bam

samtools cat -o data/Dd2comb_2D40.2_L12  data/Dd22D4_1.c4.0.2.bam data/Dd2FDK2_1.c4.0.8.bam \
    data/Dd22D4_2.c4.0.2.bam data/Dd2FDK2_2.c4.0.8.bam
samtools cat -o data/Dd2comb_2D40.4_L12  data/Dd22D4_1.c4.0.4.bam data/Dd2FDK2_1.c4.0.6.bam \
    data/Dd22D4_2.c4.0.4.bam data/Dd2FDK2_2.c4.0.6.bam
samtools cat -o data/Dd2comb_2D40.5_L12  data/Dd22D4_1.c4.0.5.bam data/Dd2FDK2_1.c4.0.5.bam \
    data/Dd22D4_2.c4.0.5.bam data/Dd2FDK2_2.c4.0.5.bam
samtools cat -o data/Dd2comb_2D40.6_L12  data/Dd22D4_1.c4.0.6.bam data/Dd2FDK2_1.c4.0.4.bam \
    data/Dd22D4_2.c4.0.6.bam data/Dd2FDK2_2.c4.0.4.bam
samtools cat -o data/Dd2comb_2D40.8_L12  data/Dd22D4_1.c4.0.8.bam data/Dd2FDK2_1.c4.0.2.bam \
    data/Dd22D4_2.c4.0.8.bam data/Dd2FDK2_2.c4.0.2.bam




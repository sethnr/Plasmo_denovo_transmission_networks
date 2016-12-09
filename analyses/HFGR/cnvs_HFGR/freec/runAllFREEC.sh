#!/bin/bash

##perl -lane 'print "bsub13 $F[0]_$F[1] runFREEC  -W 200 -B $F[4] -o $F[0]_$F[1] " ' $DISCO1/samples/Thies.samples.txt

bsub13 SM-7LV8E_1 runFREEC  -W 200 -B /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/1/SM-7LV8E/H2MGCBCXX.1.aligned.duplicates_marked.bam -o SM-7LV8E_1
bsub13 SM-7LV8E_2 runFREEC  -W 200 -B /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/2/SM-7LV8E/H2MGCBCXX.2.aligned.duplicates_marked.bam -o SM-7LV8E_2
bsub13 SM-7LV8F_1 runFREEC  -W 200 -B /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/1/SM-7LV8F/H2MGCBCXX.1.aligned.duplicates_marked.bam -o SM-7LV8F_1
bsub13 SM-7LV8F_2 runFREEC  -W 200 -B /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/2/SM-7LV8F/H2MGCBCXX.2.aligned.duplicates_marked.bam -o SM-7LV8F_2
bsub13 SM-7LV8G_1 runFREEC  -W 200 -B /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/1/SM-7LV8G/H2MGCBCXX.1.aligned.duplicates_marked.bam -o SM-7LV8G_1
bsub13 SM-7LV8G_2 runFREEC  -W 200 -B /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/2/SM-7LV8G/H2MGCBCXX.2.aligned.duplicates_marked.bam -o SM-7LV8G_2
bsub13 SM-7LV8H_1 runFREEC  -W 200 -B /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/1/SM-7LV8H/H2MGCBCXX.1.aligned.duplicates_marked.bam -o SM-7LV8H_1
bsub13 SM-7LV8H_2 runFREEC  -W 200 -B /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/2/SM-7LV8H/H2MGCBCXX.2.aligned.duplicates_marked.bam -o SM-7LV8H_2
bsub13 SM-7LV8I_1 runFREEC  -W 200 -B /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/1/SM-7LV8I/H2MGCBCXX.1.aligned.duplicates_marked.bam -o SM-7LV8I_1
bsub13 SM-7LV8I_2 runFREEC  -W 200 -B /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/2/SM-7LV8I/H2MGCBCXX.2.aligned.duplicates_marked.bam -o SM-7LV8I_2
bsub13 SM-7LV8K_1 runFREEC  -W 200 -B /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/1/SM-7LV8K/H2MGCBCXX.1.aligned.duplicates_marked.bam -o SM-7LV8K_1
bsub13 SM-7LV8K_2 runFREEC  -W 200 -B /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/2/SM-7LV8K/H2MGCBCXX.2.aligned.duplicates_marked.bam -o SM-7LV8K_2
bsub13 SM-7LV8L_1 runFREEC  -W 200 -B /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/1/SM-7LV8L/H2MGCBCXX.1.aligned.duplicates_marked.bam -o SM-7LV8L_1
bsub13 SM-7LV8L_2 runFREEC  -W 200 -B /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/2/SM-7LV8L/H2MGCBCXX.2.aligned.duplicates_marked.bam -o SM-7LV8L_2
bsub13 SM-7LV8M_1 runFREEC  -W 200 -B /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/1/SM-7LV8M/H2MGCBCXX.1.aligned.duplicates_marked.bam -o SM-7LV8M_1
bsub13 SM-7LV8M_2 runFREEC  -W 200 -B /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/2/SM-7LV8M/H2MGCBCXX.2.aligned.duplicates_marked.bam -o SM-7LV8M_2
bsub13 SM-7LV8N_1 runFREEC  -W 200 -B /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/1/SM-7LV8N/H2MGCBCXX.1.aligned.duplicates_marked.bam -o SM-7LV8N_1
bsub13 SM-7LV8N_2 runFREEC  -W 200 -B /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/2/SM-7LV8N/H2MGCBCXX.2.aligned.duplicates_marked.bam -o SM-7LV8N_2


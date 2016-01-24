#!/bin/bash

##perl -lane 'print "bsub13 $F[0]_$F[1] runFREEC -B $F[4] -o $F[0]_$F[1] " ' $DISCO1/samples/Thies.samples.txt

bsub13 SM-7LV7V_1 runFREEC -B /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/1/SM-7LV7V/H2MGCBCXX.1.aligned.duplicates_marked.bam -o SM-7LV7V_1 
bsub13 SM-7LV7V_2 runFREEC -B /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/2/SM-7LV7V/H2MGCBCXX.2.aligned.duplicates_marked.bam -o SM-7LV7V_2 
bsub13 SM-7LV7W_1 runFREEC -B /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/1/SM-7LV7W/H2MGCBCXX.1.aligned.duplicates_marked.bam -o SM-7LV7W_1 
bsub13 SM-7LV7W_2 runFREEC -B /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/2/SM-7LV7W/H2MGCBCXX.2.aligned.duplicates_marked.bam -o SM-7LV7W_2 
bsub13 SM-7LV7X_1 runFREEC -B /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/1/SM-7LV7X/H2MGCBCXX.1.aligned.duplicates_marked.bam -o SM-7LV7X_1 
bsub13 SM-7LV7X_2 runFREEC -B /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/2/SM-7LV7X/H2MGCBCXX.2.aligned.duplicates_marked.bam -o SM-7LV7X_2 
bsub13 SM-7LV7Y_1 runFREEC -B /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/1/SM-7LV7Y/H2MGCBCXX.1.aligned.duplicates_marked.bam -o SM-7LV7Y_1 
bsub13 SM-7LV7Y_2 runFREEC -B /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/2/SM-7LV7Y/H2MGCBCXX.2.aligned.duplicates_marked.bam -o SM-7LV7Y_2 
bsub13 SM-7LV7Z_1 runFREEC -B /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/1/SM-7LV7Z/H2MGCBCXX.1.aligned.duplicates_marked.bam -o SM-7LV7Z_1 
bsub13 SM-7LV7Z_2 runFREEC -B /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/2/SM-7LV7Z/H2MGCBCXX.2.aligned.duplicates_marked.bam -o SM-7LV7Z_2 
bsub13 SM-7LV81_1 runFREEC -B /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/1/SM-7LV81/H2MGCBCXX.1.aligned.duplicates_marked.bam -o SM-7LV81_1 
bsub13 SM-7LV81_2 runFREEC -B /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/2/SM-7LV81/H2MGCBCXX.2.aligned.duplicates_marked.bam -o SM-7LV81_2 
bsub13 SM-7LV82_1 runFREEC -B /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/1/SM-7LV82/H2MGCBCXX.1.aligned.duplicates_marked.bam -o SM-7LV82_1 
bsub13 SM-7LV82_2 runFREEC -B /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/2/SM-7LV82/H2MGCBCXX.2.aligned.duplicates_marked.bam -o SM-7LV82_2 
bsub13 SM-7LV83_1 runFREEC -B /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/1/SM-7LV83/H2MGCBCXX.1.aligned.duplicates_marked.bam -o SM-7LV83_1 
bsub13 SM-7LV83_2 runFREEC -B /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/2/SM-7LV83/H2MGCBCXX.2.aligned.duplicates_marked.bam -o SM-7LV83_2 
bsub13 SM-7LV84_1 runFREEC -B /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/1/SM-7LV84/H2MGCBCXX.1.aligned.duplicates_marked.bam -o SM-7LV84_1 
bsub13 SM-7LV84_2 runFREEC -B /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/2/SM-7LV84/H2MGCBCXX.2.aligned.duplicates_marked.bam -o SM-7LV84_2 
bsub13 SM-7LV85_1 runFREEC -B /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/1/SM-7LV85/H2MGCBCXX.1.aligned.duplicates_marked.bam -o SM-7LV85_1 
bsub13 SM-7LV85_2 runFREEC -B /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/2/SM-7LV85/H2MGCBCXX.2.aligned.duplicates_marked.bam -o SM-7LV85_2 
bsub13 SM-7LV86_1 runFREEC -B /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/1/SM-7LV86/H2MGCBCXX.1.aligned.duplicates_marked.bam -o SM-7LV86_1 
bsub13 SM-7LV86_2 runFREEC -B /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/2/SM-7LV86/H2MGCBCXX.2.aligned.duplicates_marked.bam -o SM-7LV86_2 
bsub13 SM-7LV87_1 runFREEC -B /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/1/SM-7LV87/H2MGCBCXX.1.aligned.duplicates_marked.bam -o SM-7LV87_1 
bsub13 SM-7LV87_2 runFREEC -B /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/2/SM-7LV87/H2MGCBCXX.2.aligned.duplicates_marked.bam -o SM-7LV87_2 
bsub13 SM-7LV88_1 runFREEC -B /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/1/SM-7LV88/H2MGCBCXX.1.aligned.duplicates_marked.bam -o SM-7LV88_1 
bsub13 SM-7LV88_2 runFREEC -B /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/2/SM-7LV88/H2MGCBCXX.2.aligned.duplicates_marked.bam -o SM-7LV88_2 
bsub13 SM-7LV89_1 runFREEC -B /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/1/SM-7LV89/H2MGCBCXX.1.aligned.duplicates_marked.bam -o SM-7LV89_1 
bsub13 SM-7LV89_2 runFREEC -B /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/2/SM-7LV89/H2MGCBCXX.2.aligned.duplicates_marked.bam -o SM-7LV89_2 
bsub13 SM-7LV8A_1 runFREEC -B /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/1/SM-7LV8A/H2MGCBCXX.1.aligned.duplicates_marked.bam -o SM-7LV8A_1 
bsub13 SM-7LV8A_2 runFREEC -B /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/2/SM-7LV8A/H2MGCBCXX.2.aligned.duplicates_marked.bam -o SM-7LV8A_2 
bsub13 SM-7LV8B_1 runFREEC -B /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/1/SM-7LV8B/H2MGCBCXX.1.aligned.duplicates_marked.bam -o SM-7LV8B_1 
bsub13 SM-7LV8B_2 runFREEC -B /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/2/SM-7LV8B/H2MGCBCXX.2.aligned.duplicates_marked.bam -o SM-7LV8B_2 
bsub13 SM-7LV8C_1 runFREEC -B /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/1/SM-7LV8C/H2MGCBCXX.1.aligned.duplicates_marked.bam -o SM-7LV8C_1 
bsub13 SM-7LV8C_2 runFREEC -B /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/2/SM-7LV8C/H2MGCBCXX.2.aligned.duplicates_marked.bam -o SM-7LV8C_2 
bsub13 SM-7LV8D_1 runFREEC -B /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/1/SM-7LV8D/H2MGCBCXX.1.aligned.duplicates_marked.bam -o SM-7LV8D_1 
bsub13 SM-7LV8D_2 runFREEC -B /seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/2/SM-7LV8D/H2MGCBCXX.2.aligned.duplicates_marked.bam -o SM-7LV8D_2 

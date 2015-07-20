#!/usr/bin/python

import sys
import argparse
from string import *
from subprocess import call, check_output
import os.path as path
from math import ceil,floor
import pysam


parser = argparse.ArgumentParser(description='SW realign variant-transformed sequence to original reference')

parser.add_argument('-v','--vcf', action="store", dest='vcf', type=str, help='vcf file', nargs='?', default=None)
parser.add_argument('-r1','-r','--ref_for_cons', action="store", dest='fasta_cons', type=str, help='fasta file (ref)', nargs='?', default=None)
parser.add_argument('-r2','-f','--fasta_for_cons', action="append", dest='fasta_align', help='fasta file (ref)', nargs='?')

parser.add_argument('-o','--out', action="store", dest='outFile', type=str, help='outFile', nargs='?', default=None)
parser.add_argument('-b','--blockSize', action="store", dest='blocksize', type=int, help='default blocksize for summary (kb)', nargs='?', default=10)
parser.add_argument('-l','--loci', action="store", dest='regions', type=str, help='regions over which to summarize', nargs='+', default=None)
args = parser.parse_args()

TMP="./tmp/"
block = args.blocksize * 1000
blocks = open(TMP+"blocks.intervals","w")

for region in args.regions:
    (chrom, locus) = region.split(":")
    (regst, regen) = locus.split("-")
    regst = int(floor(int(regst) / block)*block)
    regen = int(ceil(int(regen) / block)*block)
    for st in range(regst,regen,block):
        en = st+block
        print >>blocks,chrom+":"+str(st+1)+"-"+str(en-1)
blocks.close()

altref_command = ['java','-jar','/humgen/gsa-hpprojects/GATK/bin/GenomeAnalysisTK-3.4-0-g7e26428/GenomeAnalysisTK.jar',
                  '-T', 'FastaAlternateReferenceMaker',
                  '-R', args.fasta_cons,
                  '-o', TMP+'consensus.fasta',
#                  '-L', chrom+":"+str(st+1)+"-"+str(en),
                  '-L', TMP+'blocks.intervals',
                  '-V', args.vcf]
print >>sys.stderr, " ".join(altref_command)
check_output(altref_command)

print >>sys.stderr, args.fasta_align

quals = dict()
tags=['AS', #alignment score (no matched bases?)
      'NM', #edit distance
      'XS', #suboptimal align score
      'XF', #support from forward/reverse
      'XE'  #no supporting seeds
      ]
  
for fasta in args.fasta_align:
    bwa_command = ['bwa','bwasw',
               fasta, TMP+'consensus.fasta']
    name = path.basename(fasta).replace('.fasta','')
    samfilename = TMP+name+'.sam'
#    samfilename = TMP+"realigns.sam"
    realigns = open(samfilename,"w")
    call(bwa_command,stdout=realigns)
    realigns.close()


    #doesn't need sorting, should be in read order
    #pysam.sort("-S","-n", TMP+"realigns.sam", TMP+"realigns.sam")
    samfile = pysam.AlignmentFile(samfilename, "r")

    for a in samfile:
    #    print a.qname,
    #    print a.is_secondary,
    #    print a.qual
        q = str(a.qname)
        if (q,name,'N') not in quals:
            quals[(q,name,'N')]=0
            for tag in tags:
                quals[(q,name,tag)]=0
        quals[(q,name,'N')]+=1  # count No of aligns

        for tag,val in a.tags:
            quals[(q,name,tag)] += val

blocks = [b for b,n,t in quals]
blocks = list(set(blocks))
blocks.sort(key=int)

names = [n for b,n,t in quals]
names = list(set(names))

print '#'+args.fasta_cons
print '#'+args.vcf

print "block",
for name in names:
    print "\t".join([""] + [i+"."+name for i in ["N","AS","NM","XS"]]),
print ""
for block in blocks:
    print block,
    for name in names:
         print "\t".join(map(str,["",
                    quals[(block,name,'N')],
                    quals[(block,name,'AS')],
                    quals[(block,name,'NM')],
                    quals[(block,name,'XS')]
                    ])),
    print "" 

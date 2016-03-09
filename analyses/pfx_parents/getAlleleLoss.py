#!/usr/bin/python

import vcf
import sys
import argparse
from string import *
import collections
import copy

parser = argparse.ArgumentParser(description='get non-ref allele depth')

parser.add_argument('-v','--vcf', action="store", dest='vcfFile1', type=str, help='vcfFile1', nargs='?', default=None)

args = parser.parse_args()

vcfFile1 = open(args.vcfFile1,'r')
reader1=vcf.Reader(vcfFile1)

depthlim=20


print "\t".join(["chr","pos","type","P1","P2","P1_F1_pc","P2_F1_pc"])

P1name, P2name = reader1.samples[0:2]

print >>sys.stderr, "P1="+P1name
print >>sys.stderr, "P2="+P2name
print >>sys.stderr, "F1s="+str(reader1.samples[2:])

for rec in reader1:
    newSamples = list()
    
    nrpcs = list()
    sumDepth=0
    nrcalls=0
    F1counts = dict()
    F1called=0
    P1=rec.genotype(P1name).gt_bases
    P2=rec.genotype(P2name).gt_bases

    if P1==P2: continue

    for call in rec.samples[2:]:        
        callData= list(call.data)
        call = call.gt_bases
        if call is None: continue
        if call not in F1counts: F1counts[call]=0
        F1counts[call]+=1
        F1called+=1
    if P1 not in F1counts: P1pc='.'
    else: P1pc = round(F1counts[P1]/float(F1called),3)
    if P2 not in F1counts: P2pc='.'
    else: P2pc = round(F1counts[P2]/float(F1called),3)
    
    print >>sys.stdout, "\t".join(map(str,[
                rec.CHROM,
                rec.POS,
                rec.var_type,
                P1,P2,
                P1pc,P2pc]))

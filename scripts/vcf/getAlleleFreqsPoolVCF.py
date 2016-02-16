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

for rec in reader1:
#    calls = dict()
#    rec.samples = None
    newSamples = list()
    DP4index=rec.FORMAT.split(":").index("DP4")
    
    
    nrpcs = list()
    sumDepth=0
    nrcalls=0
    for call in rec.samples:
#        print call.data
#        print rec.FORMAT
        callData= list(call.data)
        dp4 = callData[DP4index]
        
        nrPc=-1
        callDepth=0
        if dp4 is not None:
            depths = map(int,dp4)
#            print depths
            callDepth = sum(depths)
            sumDepth += callDepth
#            print depths,(depths[2]+depths[3]), float(callDepth),
            if callDepth > 0 and callDepth > depthlim:
                nrPc = round((depths[3]+depths[2]) / float(callDepth),3)
#            else:
#                callDepth = 0
#            print nrPc
        if nrPc > 0 and callDepth > depthlim:
            nrpcs += [nrPc]
        else:
            nrpcs += [0]
    if sum(nrpcs) > 0 :
        print "\t".join(map(str,[rec.CHROM,rec.POS,rec.REF,rec.ALT[0],sumDepth] + nrpcs))

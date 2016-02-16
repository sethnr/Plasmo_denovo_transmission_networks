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


print "\t".join(["chr","pos","ref","alts","depth"]+[s+"r\t"+s+"nr" for s in reader1.samples])


for rec in reader1:
#    calls = dict()
#    rec.samples = None
    newSamples = list()
    FORMAT = rec.FORMAT.split(":")
    DP4index=-1
    ADindex=-1
    if "DP4" in FORMAT:
        DP4index=FORMAT.index("DP4")
    if "AD" in FORMAT:
        ADindex=rec.FORMAT.split(":").index("AD")
    
    nrpcs = list()
    sumDepth=0
    nrcalls=0
    for call in rec.samples:
#        print call.data
#        print rec.FORMAT
        callData= list(call.data)
        nrPc = -1
        nrC = 0
        rC = 0
        if DP4index > 0:
            dp4 = callData[DP4index]        
#            nrPc=-1
            callDepth=0
            if dp4 is not None:
                depths = map(int,dp4)
                callDepth = sum(depths)
                sumDepth += callDepth
                if callDepth > 0 and callDepth > depthlim:
                    nrPc = round((depths[3]+depths[2]) / float(callDepth),3)
                    rC = depths[0]+depths[1]
                    nrC = depths[2]+depths[3]
        elif ADindex > 0:
            ad = callData[ADindex]        
#            nrPc=-1
            callDepth=0
            if ad is not None:
                depths = map(int,ad)
                callDepth = sum(depths)
                sumDepth += callDepth
                if callDepth > 0 and callDepth > depthlim:
                    nrPc = round(depths[1] / float(callDepth),3)
                    rC = depths[0]
                    nrC = sum(depths[1:])

        if nrPc > 0 and callDepth > depthlim:
            nrpcs += [rC,nrC]
        else:
            nrpcs += [0,0]
    if sum(nrpcs) > 0 :
        print "\t".join(map(str,[rec.CHROM,rec.POS,rec.REF,"/".join(map(str,rec.ALT)),sumDepth] + nrpcs))

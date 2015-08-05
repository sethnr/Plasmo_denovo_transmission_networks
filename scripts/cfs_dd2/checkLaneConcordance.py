#!/usr/bin/python

import vcf
import sys
import argparse
from string import *

parser = argparse.ArgumentParser(description='get STR length differences from disco VCF')

parser.add_argument('-v','--vcf', action="store", dest='vcfFile1', type=str, help='vcfFile1', nargs='?', default=None)

parser.add_argument('-o','--out', action="store", dest='outFile', type=str, help='make combined outfile named <outfile> instead of separate files', nargs='?', default=None)
#parser.add_argument('-i','--defaultSize', action="store", dest='defSize', type=int, help='default size for STRs', nargs='?', default=0)


args = parser.parse_args()

vcfFile1 = open(args.vcfFile1,'r')
reader1=vcf.Reader(vcfFile1)

vcfoutF1 = replace(args.vcfFile1,'.vcf','.lanes.vcf')
vcfoutF1 = replace(vcfoutF1,'.vcf.gz','.lanes.vcf')
print >>sys.stdout, "VCFOUT1",vcfoutF1
vcfoutF1 = open(vcfoutF1,'w')
vcfout1=vcf.Writer(vcfoutF1,reader1)

samples = set()
for s in reader1.samples:
    sample, lane = s.split('_')
    samples.add(sample)
samples = list(samples)
samples = sorted(samples)

print "\t".join(["chr","pos","alleles","type"]+samples)
for rec in reader1:
    calls = dict()
    for call in rec.samples:
        sample, lane = call.sample.split('_')
        calls[(sample,int(lane))] = call.gt_type
#    samples = list(set([s for s,l in calls]))
    for s in samples:
        if (s,1) not in calls: calls[(s,1)] = None
        if (s,2) not in calls: calls[(s,2)] = None

        if calls[(s,1)] is None and calls[(s,2)] is None:
            calls[(s,-1)] = "."
        elif calls[(s,1)] is None:
            calls[(s,-1)] = str(calls[(s,2)])+"*"
        elif calls[(s,2)] is None:        
            calls[(s,-1)] = str(calls[(s,1)])+"*"
        elif calls[(s,1)] == calls[(s,2)]:
            calls[(s,-1)] = str(calls[(s,1)])
        else:
            calls[(s,-1)] = "!!"
    if rec.is_indel:
        print "\t".join([rec.CHROM,str(rec.POS),str(len(rec.alleles)),rec.var_subtype]+
                        [calls[(s,-1)] for s in samples])



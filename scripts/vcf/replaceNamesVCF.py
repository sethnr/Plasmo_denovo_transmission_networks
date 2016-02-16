#!/usr/bin/python

import vcf
import sys
import argparse
from string import *
import collections
import copy

parser = argparse.ArgumentParser(description='get STR length differences from disco VCF')

parser.add_argument('-v','--vcf', action="store", dest='vcfFile1', type=str, help='vcfFile1', nargs='?', default=None)
parser.add_argument('-f','--lookup', action="store", dest='lookup', type=str, help='sample name lookup file', nargs='?', default=None)

args = parser.parse_args()

vcfFile1 = open(args.vcfFile1,'r')
reader1=vcf.Reader(vcfFile1)

vcfoutF1 = replace(args.vcfFile1,'.vcf','.RENAME.vcf')
vcfoutF1 = replace(vcfoutF1,'.vcf.gz','.vcf')
print >>sys.stdout, "VCFOUT1",vcfoutF1
vcfoutF1 = open(vcfoutF1,'w')

vcfsamples = reader1.samples
samples = list()
discards = list()
sampleKey = dict()
if args.lookup is not None:
    lfile = open(args.lookup,'r')
    for l in lfile:
        l = l.rstrip()
        lookup = l.split("\t")
#        subsamp, sample = lookup[0],lookup[1]
        subsamp = lookup[0]
        sample = lookup[1]

        print subsamp
        print sample
        sampleKey[subsamp]=sample
    pass
else:
    print >>sys.stderr, "new sample names not provided"
    exit(1)

for oldsamp in vcfsamples:
    if oldsamp not in sampleKey:
        sampleKey[oldsamp] = oldsamp
#    print >>sys.stderr, "  "+oldsamp+" - "+sampleKey[oldsamp],
    newsamp = sampleKey[oldsamp]
    if newsamp not in samples:
        samples += [newsamp]
    else:
        discards += [oldsamp]
    
for key in vcfsamples:
    print >>sys.stderr, "  "+key+" - "+sampleKey[key],
#    if key not in vcfsamples: 
#        print >>sys.stderr, "\tnot in vcf"
    if key in discards:
        print >>sys.stderr, "\tduplicate"
    else:
        print >>sys.stderr,""
        
reader2 = copy.copy(reader1)
reader2.samples = samples
print len(reader1.samples)
print len(reader2.samples)

vcfout1=vcf.Writer(vcfoutF1,reader2)



for rec in reader1:
#    calls = dict()
#    rec.samples = None
    newSamples = list()
    newSampleNames = list()
    for call in rec.samples:
        #print >> sys.stderr,  call.sample
        if sampleKey[call.sample] not in newSampleNames:
            newSampleNames += [sampleKey[call.sample]]
            call.sample=sampleKey[call.sample]
            newSamples += [call]
    #print >>sys.stderr, len(rec.samples), len(newSamples)
    rec.samples = newSamples
    vcfout1.write_record(rec)

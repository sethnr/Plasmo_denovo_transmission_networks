#!/usr/bin/python

import vcf
import sys
import argparse
from string import *

parser = argparse.ArgumentParser(description='get STR length differences from disco VCF')

parser.add_argument('-v','--vcf', action="store", dest='vcfFile', type=str, help='vcfFile', nargs='?', default=None)
parser.add_argument('-i','--defaultSize', action="store", dest='defSize', type=int, help='default size for STRs', nargs='?', default=0)


args = parser.parse_args()
print >>sys.stderr, args.vcfFile

vcfFile = open(args.vcfFile,'r')
reader=vcf.Reader(vcfFile)

print "\t".join(["#CHROM","POS"]),
print "\t".join(reader.samples)
    
for record in reader:
    refAllele = record.alleles[0]
    STRlens = dict()
    
    for allele in record.alleles:
        STRlen = len(allele) - len(refAllele) + args.defSize
        STRlens[str(allele)] = str(STRlen)
#    print record.alleles
#    print STRlens
    print "\t".join([str(record.CHROM),str(record.POS)]),
    for call in record.samples:
        callSTRlens = []
        if call.gt_type is None:
            print "\t.",
        else:
            print "\t",
            for allele in split(call.gt_bases,"/"):
                callSTRlens += [STRlens[allele]]
            print "/".join(callSTRlens),
    print ""

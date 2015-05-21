#!/usr/bin/python

import vcf
import sys
import argparse
from string import *

parser = argparse.ArgumentParser(description='get STR length differences from disco VCF')

parser.add_argument('-v','--vcf', action="store", dest='vcfFile', type=str, help='vcfFile', nargs='?', default=None)
parser.add_argument('-o','--out', action="store", dest='outFile', type=str, help='outFile', nargs='?', default=None)
parser.add_argument('-i','--defaultSize', action="store", dest='defSize', type=int, help='default size for STRs', nargs='?', default=0)


args = parser.parse_args()
print >>sys.stderr, args.vcfFile


outfile = sys.stderr
if args.outFile is not None:
    outfile = open(args.outFile,'w')

    
vcfFile = open(args.vcfFile,'r')
reader=vcf.Reader(vcfFile)

print >>outfile, "\t".join(["#CHROM","POS"]),
print >>outfile, "\t".join(reader.samples)
    
for record in reader:
    refAllele = record.alleles[0]
    STRlens = dict()
    
    for allele in record.alleles:
        STRlen = len(allele) - len(refAllele) + args.defSize
        STRlens[str(allele)] = str(STRlen)
#    print record.alleles
#    print STRlens

    print >>outfile, "\t".join([str(record.CHROM),str(record.POS)]),
    for call in record.samples:
        callSTRlens = []
        if call.gt_type is None:
            print >>outfile, "\t.",
        else:
            print >>outfile, "\t",
            for allele in split(call.gt_bases,"/"):
                callSTRlens += [STRlens[allele]]
            print >>outfile, "/".join(callSTRlens),
    print >>outfile, ""

sys.exit(0)

#!/usr/bin/python

import vcf
import sys
import argparse
from string import *
import collections
import copy

parser = argparse.ArgumentParser(description='get allele numbers table')

parser.add_argument('-v','--vcf', action="store", dest='vcfFile1', type=str, help='vcfFile1', nargs='?', default=None)

args = parser.parse_args()

vcfFile1 = open(args.vcfFile1,'r')
reader1=vcf.Reader(vcfFile1)

print "\t".join(["chr","pos","type","alleles"]+reader1.samples)
for rec in reader1:
    if rec.num_unknown > 2: continue

    alleles = "/".join([rec.REF] + map(str,rec.ALT))
    vartype="SNP";
    if rec.is_indel: vartype="INDEL"
    print "\t".join(map(str,[rec.CHROM,rec.POS,vartype,alleles])),
    GTindex=rec.FORMAT.split(":").index("GT")
#    print GTindex
    for call in rec.samples:
        callData=call.data[GTindex]
        if callData is not None:
            print "\t"+callData,
        else:
            print "\t.",
    print ""

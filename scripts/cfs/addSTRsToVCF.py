#!/usr/bin/python

import vcf
import sys
import argparse
from string import *
from os import path
import re
import bz2
from math import floor

parser = argparse.ArgumentParser(description='take data from mreps output file, add to VCF')

parser.add_argument('-v','--vcf', action="store", dest='vcfFile', type=str, help='vcfFile', nargs='?', default=None)
parser.add_argument('-s','--str', action="store", dest='strFile', type=str, help='outFile', nargs='?', default=None)
args = parser.parse_args()


## get complexity of string (cf to compressed size)
def stringComplex(mystring):
    comp = bz2.compress(mystring)
    bcomplex = float(len(mystring))/len(comp)
    return bcomplex



print >>sys.stderr, "reading in STR file "+args.strFile
strfile = open(args.strFile,'r')
STR = dict()
for line in strfile:
    F = line.split()
    if len(F) > 8:
        (chrom, start, arrow, end, colon, size, period, exponent, error) = F[:9]
        components = F[9:]
        seq = "".join(components)
        complexity = round(stringComplex(seq),3)
        period = int(re.sub('\D','',period))
        exponent = float(re.sub('[\[\]]','',exponent))

        fullreps = int(floor(exponent))
        fullrepseq = "".join(components[:fullreps])
        Apc = sum([1 for i in fullrepseq if i in ('a','A')])/float(len(fullrepseq))
        Tpc = sum([1 for i in fullrepseq if i in ('t','T')])/float(len(fullrepseq))
        Cpc = sum([1 for i in fullrepseq if i in ('c','C')])/float(len(fullrepseq))
        Gpc = sum([1 for i in fullrepseq if i in ('g','G')])/float(len(fullrepseq))
        basePC = ":".join([str(round(i,2)) for i in [Apc,Tpc,Cpc,Gpc]])

        start = int(start)
        end = int(end)
        size = int(size)
        
        for p in range(start,end):
            if (chrom,p) in STR:
                (start2, end2, size2, period2, exponent2,complexity,basePC) = STR[(chrom,p)]
                if size > size2:
                    STR[(chrom,p)] = (start, end, size, period, exponent,complexity,basePC)
            else:
                STR[(chrom,p)] = (start, end, size, period, exponent,complexity,basePC)
                    
vcfFile = open(args.vcfFile,'r')
reader=vcf.Reader(vcfFile)

vcfin = path.basename(args.vcfFile)
vcfout = replace(vcfin,'.vcf','.STRs.vcf')
vcfout = replace(vcfout,'.vcf.gz','.STRs.vcf')

vcfout = open(vcfout,'w')
writer=vcf.Writer(vcfout,reader)

print >>sys.stderr, vcfFile
for rec in reader:
    if (rec.CHROM,rec.POS) in STR:
        (start, end, size, period, exponent,complexity,basePC) = STR[(rec.CHROM,rec.POS)]
        rec.INFO['STR']=rec.CHROM+":"+str(start)+"-"+str(end)
        rec.INFO['STRP']=period
        rec.INFO['STRE']=exponent
        rec.INFO['STRS']=size
        rec.INFO['STRC']=complexity
        rec.INFO['STRatcg']=basePC
    vartype="SNP"
    if rec.is_indel: vartype="INDEL"
    rec.INFO['type']=vartype
        
    writer.write_record(rec)


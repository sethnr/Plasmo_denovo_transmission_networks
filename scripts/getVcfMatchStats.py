#!/usr/bin/python

import vcf
import sys
import argparse
from string import *
import bz2
from os import path

parser = argparse.ArgumentParser(description='get STR length differences from disco VCF')

parser.add_argument('-v','--vcf', action="store", dest='vcfFile1', type=str, help='vcfFile1', nargs='?', default=None)

args = parser.parse_args()
vcfFile1 = open(args.vcfFile1,'r')
reader=vcf.Reader(vcfFile1)

vcfname = path.basename(args.vcfFile1)
vcfname = path.splitext(vcfname)[0]
for rec in reader:
    alleles = [str(i) for i in rec.alleles]
    varlen=1
    bcomplex = 1
    if rec.is_indel:        
        maxallele = max(alleles, key=len)
        varlen = len(maxallele)
        comp = bz2.compress(maxallele)
#        print len(maxallele), maxallele
#        print len(comp)
        bcomplex = float(len(maxallele))/len(comp)
    noAlls = str(len(alleles))
    bcomplex = str(bcomplex)
    varlen = str(varlen)
    type = rec.INFO['TYPE'][0]
    con = rec.INFO['CON'][0]
    qual = str(rec.QUAL)

    STR=""
    STRp='0'
    STRe='0'
    STRs='0'
    if 'STR' in rec.INFO:
        STR=rec.INFO['STR'][0]
        STRp=rec.INFO['STRP'][0]
        STRe=rec.INFO['STRE'][0]
        STRs=rec.INFO['STRS'][0]  
    
    print >>sys.stdout,"\t".join((vcfname,
                                  rec.CHROM,
                                  str(rec.POS),
                                  type,
                                  con,
                                  qual,
                                  varlen,
                                  noAlls,
                                  bcomplex,
                                  STR,
                                  STRp,
                                  STRe,
                                  STRs
                                  ))

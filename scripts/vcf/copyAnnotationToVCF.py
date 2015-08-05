#!/usr/bin/python

import vcf
import sys
import argparse
from string import *

parser = argparse.ArgumentParser(description='from CHR : POS : 1/0 file, add true/false to VCF')

parser.add_argument('-v1','-vcfFrom', action="store", dest='vcfFrom', type=str, help='vcf file', nargs='?', default=None)
parser.add_argument('-v2','-vcfTo', action="store", dest='vcfTo', type=str, help='valueFile', nargs='?', default=None)

parser.add_argument('-i','--infoname', action="store", dest='valname', type=str, help='make outfile named <vcfTo>.<infoname>.vcf. ', nargs='?', default="VALUE")
#parser.add_argument('-i','--defaultSize', action="store", dest='defSize', type=int, help='default size for STRs', nargs='?', default=0)


args = parser.parse_args()

vcfFileFrom = open(args.vcfFrom,'r')
vcfFrom=vcf.Reader(vcfFileFrom)

vcfFileTo = open(args.vcfTo,'r')
vcfTo=vcf.Reader(vcfFileTo)

vcfoutF = replace(args.vcfTo,'.vcf','.'+args.valname+'.vcf')
vcfoutF = replace(vcfoutF,'.vcf.gz','.'+args.valname+'.vcf')
#print >>sys.stdout, "VCFOUT1",vcfoutF1
vcfoutF = open(vcfoutF,'w')
vcfout=vcf.Writer(vcfoutF,vcfTo)

vals = dict()
for rec in vcfFrom:
    #if p == 451480: print c,p,v,bool(v) 
    vals[(rec.CHROM,rec.POS)] = rec.INFO[args.valname]

for rec in vcfTo:
    #if rec.is_indel:
    if (rec.CHROM,rec.POS) in vals:
        rec.INFO[args.valname]=vals[(rec.CHROM,rec.POS)]
    vcfout.write_record(rec)


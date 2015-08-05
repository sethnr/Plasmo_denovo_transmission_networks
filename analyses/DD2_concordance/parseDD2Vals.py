#!/usr/bin/python

import sys
import argparse
from string import *
from subprocess import call, check_output
import os.path as path
from math import floor,ceil


import numpy as np
import vcf
#from pysam import TabixFile
import tabix

parser = argparse.ArgumentParser(description='get STR length differences from disco VCF')

parser.add_argument('-v','--vcf', action="store", dest='vcfFile', type=str, help='vcf file', nargs='?', default=None)
parser.add_argument('-i','--infoname', action="store", dest='valname', type=str, help='make outfile named <vcfTo>.<infoname>.vcf. ', nargs='+', default=None)

parser.add_argument('-b','--blockSize', action="store", dest='blocksize', type=int, help='default blocksize for summary (kb)', nargs='?', default=1)

parser.add_argument('-r','--region', action="store", dest='regions', type=str, help='regions over which to summarize', nargs='+', default=None)

args = parser.parse_args()
outfile = sys.stdout

if not path.isfile(args.vcfFile+".tbi"):
    #make tabix index for vcf if not tabixed
    if path.splitext(args.vcfFile) != '.gz':
        call(["bgzip",args.vcfFile])
        args.vcfFile += '.gz'
    call(["tabix","-p","vcf",args.vcfFile])
vcfFile = open(args.vcfFile,'r')
reader=vcf.Reader(vcfFile)


def _getVal(rec,valname):
    retval = ""
    if valname in rec.INFO:
        val = rec.INFO[valname]
        if type(val)is bool:
            retval = valname
        elif type(val) is str:
            retval = rec.INFO[valname]
        elif type(val) is list:
            retval = ",".join(rec.INFO[valname])
    return retval

infovals = ['CALLABLE','MATCHED','DD2CONC','DD2CONS','DD2ALL','STR','STRatcg']

print "\t".join(["CHROM","POS","ALLELES","DD2CALLRATE"]+infovals)

for rec in reader:
    vartype = "UNKNOWN"
    if rec.is_indel: vartype="INDEL"
    elif rec.is_snp: vartype="SNP"
    
    vals = []
    for valname in infovals:
        vals += [_getVal(rec,valname)]

    calls = 0
    counts = 0
    for call in rec.samples[:8]:
        counts +=1
        if call.gt_type is not None:
            calls +=1
    call_rate = float(calls)/counts

    print "\t".join(map(str,[
                rec.CHROM,rec.POS,vartype,
                len(rec.ALT),
                call_rate]
#                _getVal(rec,'CALLABLE'),
#                _getVal(rec,'MATCHED'),
#                _getVal(rec,'DD2CONC'),
#                _getVal(rec,'DD2CONS'),
#                _getVal(rec,'DD2ALL'),
#                _getVal(rec,'STR'),
#                _getVal(rec,'STRatcg')
                + vals
                ))

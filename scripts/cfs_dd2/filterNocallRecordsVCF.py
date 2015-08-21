#!/usr/bin/python

import vcf
import sys
import argparse
from string import *
import collections
import numpy as np

parser = argparse.ArgumentParser(description='get STR length differences from disco VCF')

parser.add_argument('-v','--vcf', action="store", dest='vcfFile1', type=str, help='vcfFile1', nargs='?', default=None)
parser.add_argument('-c','--callRate', action="store", dest='callRate', type=float, help='minimum call rate - 1  all genotypes must be called', nargs='?', default=1)


args = parser.parse_args()

vcfFile1 = open(args.vcfFile1,'r')
reader1=vcf.Reader(vcfFile1)

_Info = collections.namedtuple('Info', ['id', 'num', 'type', 'desc', 'source', 'version'])
_Filter = collections.namedtuple('Filter', ['id', 'desc'])
_Alt = collections.namedtuple('Alt', ['id', 'desc'])
_Format = collections.namedtuple('Format', ['id', 'num', 'type', 'desc'])

reader1.infos['CR'] = _Info(id='CR', num=1, type='Float', desc='call-rate across all samples', source="filterNocallRecords", version=1)

fname = 'NoCall'+str(1-args.callRate)
reader1.filters[fname] = _Filter(id=fname, desc='no-call rate greater than '+str(1-args.callRate))


vcfoutF1 = replace(args.vcfFile1,'.vcf','.'+fname+'.vcf')
vcfoutF1 = replace(vcfoutF1,'.vcf.gz','.'+fname+'.vcf')
vcfoutF1 = open(vcfoutF1,'w')
vcfout1=vcf.Writer(vcfoutF1,reader1)

samples = reader1.samples

#maxAlleles=len(samples)
maxAlleles=50


makeCallData = vcf.model.make_calldata_tuple(("GT","ALTP","REFP","GP"))


for rec in reader1:
    calls = dict()
    altCount = len(rec.ALT)
    if rec.call_rate < args.callRate:
        if len(rec.FILTER) == -1:
            rec.FILTER[1]=fname
        else:
            rec.FILTER += [fname]
            
        
    
    rec.INFO['CR'] = args.callRate
    vcfout1.write_record(rec)


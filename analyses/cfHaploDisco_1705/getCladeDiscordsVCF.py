#!/usr/bin/python

import vcf
import sys
import argparse
from string import *
import collections
parser = argparse.ArgumentParser(description='get STR length differences from disco VCF')

parser.add_argument('-v','--vcf', action="store", dest='vcfFile1', type=str, help='vcfFile1', nargs='?', default=None)
parser.add_argument('-s','--samples', action="append", dest='samples', type=str, help='samples to calc across', nargs='?', default=None)
parser.add_argument('-p','--prefixVars', action="store", dest='prefix', type=str, help='prefix for outfile / info / filters', nargs='?', default="")
parser.add_argument('-F','--noFilter', action="store_false", dest='filter', help='do not add filter if MAF > 1', default=True)
parser.add_argument('-M','--ignoreMinors', action="store_true", dest='ignoreMinors', help='ignore minor allele calls when looking at concordance', default=False)

parser.add_argument('-o','--out', action="store", dest='outFile', type=str, help='make combined outfile named <outfile> instead of separate files', nargs='?', default=None)
#parser.add_argument('-i','--defaultSize', action="store", dest='defSize', type=int, help='default size for STRs', nargs='?', default=0)


args = parser.parse_args()

vcfFile1 = open(args.vcfFile1,'r')
reader1=vcf.Reader(vcfFile1)

_Info = collections.namedtuple('Info', ['id', 'num', 'type', 'desc', 'source', 'version'])
_Filter = collections.namedtuple('Filter', ['id', 'desc'])
_Alt = collections.namedtuple('Alt', ['id', 'desc'])
_Format = collections.namedtuple('Format', ['id', 'num', 'type', 'desc'])

reader1.infos[args.prefix+'MAF'] = _Info(id=args.prefix+'MAF', num=1, type='Float', desc=args.prefix+' missingness', source="checkSampleConsistency.py", version=1)
reader1.infos[args.prefix+'MISS'] = _Info(id=args.prefix+'MISS', num=1, type='Float', desc=args.prefix+' missingness', source="checkSampleConsistency.py", version=1)
if args.filter:
    reader1.filters[args.prefix+'Incons'] = _Filter(id=args.prefix+'Incons', desc=args.prefix+' samples not consistent')


vcfoutF1 = replace(args.vcfFile1, '.vcf.gz',  '.vcf')
vcfoutF1 = replace(vcfoutF1,            '.vcf',       '.'+args.prefix+'.vcf')
#print >>sys.stdout, "VCFOUT1",vcfoutF1
vcfoutF1 = open(vcfoutF1,'w')
vcfout1=vcf.Writer(vcfoutF1,reader1)

if args.samples is None:    
    samples = set()
    for s in reader1.samples:
        sample, lane = s.split('_')
        samples.add(sample)
    samples = list(samples)
    samples = sorted(samples)
else:
    samples = args.samples

print samples


for rec in reader1:
    calls = dict()

    #calculate call concordance
    conc=None
    var = rec.CHROM+":"+str(rec.POS)
            
    #setup alleles count hash
    alleles = dict()
    alleles[rec.REF]=0
    for alt in rec.ALT:
        alleles[str(alt)]=0

    genos = 0
    missing = 0
    for sample in samples:
#        for lane in ('1','2'):
#        call = rec.genotype(sample+"_"+lane)
        call = rec.genotype(sample)
        gt = call.gt_bases
        if gt is not None:
            if gt.find('/')>=0:
                a1, a2 = gt.split('/')
                alleles[a1]+=1
                genos+=1
                if not args.ignoreMinors:
                    alleles[a2]+=1
                    genos +=1
            else:
                a1 = gt
                alleles[a1]+=1
                genos+=1
              
        else:
            missing +=1
    major = ''
    majorC = 0
    acopy = [a for a in alleles]
    for a in acopy:
        if alleles[a] == 0: del alleles[a] 
        else:
            if alleles[a] > majorC:
                majorC = alleles[a]
                major=a
    MAF = -1
    if genos > 0:
        MAF = (genos-majorC)/float(genos)
    MISS = missing/float(len(samples)*2)
    rec.INFO[args.prefix+'MAF'] = MAF
    rec.INFO[args.prefix+'MISS'] = MISS
    
    if args.filter and MAF > 0:
        if len(rec.FILTER) == -1:
            rec.FILTER[1]=args.prefix+'Incons'
        else:
            rec.FILTER += [args.prefix+'Incons']

#    if rec.is_snp: rec.INFO['type']='SNP'
#    else: rec.INFO['type']='INDEL'
    vcfout1.write_record(rec)


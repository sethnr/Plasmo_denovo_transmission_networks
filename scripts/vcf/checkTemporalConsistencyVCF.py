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


vcfoutF1 = replace(args.vcfFile1,'.vcf','.'+args.prefix+'.vcf')
vcfoutF1 = replace(vcfoutF1,'.vcf.gz','.'+args.prefix+'.vcf')
#print >>sys.stdout, "VCFOUT1",vcfoutF1
#vcfoutF1 = open(vcfoutF1,'w')
#vcfout1=vcf.Writer(vcfoutF1,reader1)

if args.samples is None:    
    samples = list()
    for s in reader1.samples:
        if s.find('_') > -1:
            sample, lane = s.split('_')
        else:
            sample = s
        samples.add(sample)
    samples = list(samples)
# Don't sort - assume given in temporal order
#    samples = sorted(samples)
else:
    samples = args.samples

print samples


for rec in reader1:
    calls = dict()

    #calculate call concordance
    conc=None
    var = rec.CHROM+":"+str(rec.POS)
            
    #setup alleles count hash
    alleles = list()
    alCount=dict()
    alleles += [rec.REF]
    alCount[0]=0
    i=0
    for alt in rec.ALT:
        i+=1
        alleles+=[str(alt)]
        alCount[i]=0

    genos = 0
    missing = 0
    muts=0
    startGeno = None
    currentGeno = ''
    status="none"

    GTindex=rec.FORMAT.split(":").index("GT")

    for sample in samples:
        call = rec.genotype(sample)
        gt = call.data[GTindex]
        if gt is not None:
            if gt.find('/') > -1:
                a1, a2 = gt.split('/')
                alCount[int(a1)]+=1
                alCount[int(a2)]+=1              
            else:
                a1 = gt
                alCount[int(a1)]+=1
            if startGeno is None:
                startGeno = gt
#                currentGeno = gt
            elif gt != currentGeno:
                muts += 1
                if gt != startGeno:
                    status="mut"+str(muts)
                else:
                    status="revertant"
            currentGeno = gt
#           if not args.ignoreMinors:
#              alleles[a2]+=1
            genos +=1
        else:
            missing +=1
    
    #print debug to screen
    vartype="SNP"
    if rec.is_indel: vartype="INDEL"
    print >>sys.stdout, "\t".join(map(str,
                                      [rec.CHROM,rec.POS,vartype,
                                       alleles])),
    for sample in samples:
        call = rec.genotype(sample)
        gt = call.data[GTindex]
        if gt is None: gt = '.'
        print >>sys.stdout,"\t"+str(gt),
    if startGeno is not None:
        allGenos = sum(alCount.values())
        mutCount = float(allGenos-alCount[0])/allGenos
#        mutCount = len(samples)-missing-alCount[int(startGeno)]
    else: mutCount = 0
    print status+"\t"+str(round(mutCount,3))
#    print status
        

#    if rec.is_snp: rec.INFO['type']='SNP'
#    else: rec.INFO['type']='INDEL'
#    vcfout1.write_record(rec)


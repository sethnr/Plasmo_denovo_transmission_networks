#!/usr/bin/python

import vcf
import sys
import argparse
from string import *
import collections
parser = argparse.ArgumentParser(description='check for concordance with single-sample truth set (nucmer ref comparisons?)')

parser.add_argument('-v1','--vreal', action="store", dest='vcfFile1', type=str, help='vcfFile1: real DD2', nargs='?', default=None)
parser.add_argument('-v2','--vfake', action="store", dest='vcfFile2', type=str, help='vcfFile2: fake DD2', nargs='?', default=None)

parser.add_argument('-s','--samples', action="append", dest='samples', type=str, help='samples to calc across', nargs='?', default=None)
parser.add_argument('-p','--prefixVars', action="store", dest='prefix', type=str, help='prefix for outfile / info / filters', nargs='?', default="DD2")
parser.add_argument('-F','--noFilter', action="store_false", dest='filter', help='do not add filter if MAF > 1', default=True)
parser.add_argument('-M','--ignoreMinors', action="store_true", dest='ignoreMinors', help='ignore minor allele calls when looking at concordance', default=False)

args = parser.parse_args()

vcfFile1 = open(args.vcfFile1,'r')
realReader=vcf.Reader(vcfFile1)

vcfFile2 = open(args.vcfFile2,'r')
fakeReader=vcf.Reader(vcfFile2)

_Info = collections.namedtuple('Info', ['id', 'num', 'type', 'desc', 'source', 'version'])
_Filter = collections.namedtuple('Filter', ['id', 'desc'])
_Alt = collections.namedtuple('Alt', ['id', 'desc'])
_Format = collections.namedtuple('Format', ['id', 'num', 'type', 'desc'])

realReader.infos[args.prefix+'CONC'] = _Info(id=args.prefix+'CONC', num=1, type='Float', desc=args.prefix+' concordance with fake sample', source="checkSampleConcordance.py::"+args.vcfFile2, version=1)
if args.filter:
    realReader.filters[args.prefix+'Discord'] = _Filter(id=args.prefix+'Discord', desc=args.prefix+' found in DD2 ref, but allele not concordant')
    realReader.filters[args.prefix+'RefAbsent'] = _Filter(id=args.prefix+'RefAbsent', desc=args.prefix+' var not found in DD2 ref')


vcfoutF1 = replace(args.vcfFile1,'.vcf','.'+args.prefix+'CONC.vcf')
vcfoutF1 = replace(vcfoutF1,'.vcf.gz','.'+args.prefix+'CONC.vcf')
#print >>sys.stdout, "VCFOUT1",vcfoutF1

vcfoutF1 = open(vcfoutF1,'w')
vcfout1=vcf.Writer(vcfoutF1,realReader)

samples = set()
for s in realReader.samples:
    sample, lane = s.split('_')
    samples.add(sample)
samples = list(samples)
samples = sorted(samples)

fakecall = dict()
for rec in fakeReader:
    #if rec.is_indel:
    var = rec.CHROM+":"+str(rec.POS)
    fakecall[var] = rec.REF+"/"+str(rec.ALT[0])
#    print var+" "+rec.REF+"/"+str(rec.ALT[0])

#DD2s=["SM-7LV8O","SM-7LV8P","SM-7LV8Q","SM-7LV8R"]

#print "\t".join(["chr","pos","alleles","type","conc","cons","all"]+samples)

fCounts=dict()

for rec in realReader:
    calls = dict()

    #calculate call concordance
    conc=None
    var = rec.CHROM+":"+str(rec.POS)
    
    alleles = dict()
    alleles[rec.REF]=0
    for alt in rec.ALT:
        alleles[str(alt)]=0


#    for call in rec.samples:
#        sample, lane = call.sample.split('_')
#        calls[(sample,int(lane))] = call.gt_bases
    genos =0
    missing = 0
    
    for sample in samples:
        for lane in ('1','2'):
            call = rec.genotype(sample+"_"+lane)
            gt = call.gt_bases
            if gt is not None:
                a1, a2 = gt.split('/')
                alleles[a1]+=1
                genos+=1
                alleles[a2]+=1
                genos +=1
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

    conc = False
    #if has any alt alleles (some can be refcalls only)
#    print alleles
    if bool(alleles):
        valkey = [(alleles[k],k) for k in alleles]
        valkey.sort()
#        print rec.POS,alleles, valkey
        sortkeys = [k for (v,k) in valkey]
        if rec.REF in sortkeys:
            sortkeys.remove(rec.REF) #deal only with alt calls
    
        #make call from ref+most common alt
        if len(sortkeys) > 0:
            realCall = rec.REF+"/"+sortkeys[0]
        else:
            realCall = rec.REF+"/"+rec.REF
        if var in fakecall:
            if realCall == fakecall[var]:
                conc = True
        else:
            conc=None
            


    if rec.is_snp: 
        rec.INFO['type']='SNP'
    else: 
        rec.INFO['type']='INDEL'

    #GET COUNTS OF WHICH FILTERS REMOVE DISCORDANT ALLELES
    if len(rec.FILTER) == 0:
        if ('PASS',rec.INFO['type'],conc) not in fCounts: 
            fCounts[('PASS',rec.INFO['type'],conc)] = 0
        fCounts[('PASS',rec.INFO['type'],conc)] += 1
   
    else:
        for f in rec.FILTER:
            if (f,rec.INFO['type'],conc) not in fCounts: 
                fCounts[(f,rec.INFO['type'],conc)] = 0
            fCounts[(f,rec.INFO['type'],conc)] += 1
           
    if args.filter and conc is None:
        if len(rec.FILTER) == -1:
            rec.FILTER[1]=args.prefix+'refAbsent'
        else:
            rec.FILTER += [args.prefix+'RefAbsent']
    elif args.filter and not conc:
        if len(rec.FILTER) == -1:
            rec.FILTER[1]=args.prefix+'Discord'
        else:
            rec.FILTER += [args.prefix+'Discord'] 
        
    vcfout1.write_record(rec)

filters = [f for f,t,c in fCounts]
filters = list(set(filters))
filters.sort()
print "\t".join(["FILTER","S_CONC","S_DISC","S_MISS","I_CONC","I_DISC","I_MISS"])
for f in filters:
    if (f,'SNP',True) not in fCounts: fCounts[(f,'SNP',True)]=0
    if (f,'SNP',False) not in fCounts: fCounts[(f,'SNP',False)]=0
    if (f,'INDEL',True) not in fCounts: fCounts[(f,'INDEL',True)]=0
    if (f,'INDEL',False) not in fCounts: fCounts[(f,'INDEL',False)]=0

    print "\t".join(map(str,[f, fCounts[(f,'SNP',True)],fCounts[(f,'SNP',False)],fCounts[(f,'SNP',None)],
                    fCounts[(f,'INDEL',True)],fCounts[(f,'INDEL',False)],fCounts[(f,'INDEL',None)]]))
                    
    
    

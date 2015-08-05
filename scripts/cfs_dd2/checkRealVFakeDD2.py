#!/usr/bin/python

import vcf
import sys
import argparse
from string import *

parser = argparse.ArgumentParser(description='get STR length differences from disco VCF')

parser.add_argument('-v1','--vreal', action="store", dest='vcfFile1', type=str, help='vcfFile1: real DD2', nargs='?', default=None)
parser.add_argument('-v2','--vfake', action="store", dest='vcfFile2', type=str, help='vcfFile2: fake DD2', nargs='?', default=None)

parser.add_argument('-o','--out', action="store", dest='outFile', type=str, help='make combined outfile named <outfile> instead of separate files', nargs='?', default=None)
#parser.add_argument('-i','--defaultSize', action="store", dest='defSize', type=int, help='default size for STRs', nargs='?', default=0)


args = parser.parse_args()

vcfFile1 = open(args.vcfFile1,'r')
realReader=vcf.Reader(vcfFile1)

vcfFile2 = open(args.vcfFile2,'r')
fakeReader=vcf.Reader(vcfFile2)

vcfoutF1 = replace(args.vcfFile1,'.vcf','.DD2c.vcf')
vcfoutF1 = replace(vcfoutF1,'.vcf.gz','.DD2c.vcf')
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

DD2s=["SM-7LV8O","SM-7LV8P","SM-7LV8Q","SM-7LV8R"]

print "\t".join(["chr","pos","alleles","type","conc","cons","all"]+samples)
for rec in realReader:
    calls = dict()

    #calculate call concordance
    conc=None
    var = rec.CHROM+":"+str(rec.POS)
            
    if var in fakecall:
        conc=False
        for alt in rec.ALT:
            alleles = rec.REF+"/"+str(alt)
#            print "THIS",var, fakecall[var], alleles
            if fakecall[var] == alleles:
                conc=True
    else: 
        pass
#        print rec.CHROM,rec.POS,"!"

    for call in rec.samples:
        sample, lane = call.sample.split('_')
        calls[(sample,int(lane))] = call.gt_type

#    samples = list(set([s for s,l in calls]))
    for s in samples:
        if (s,1) not in calls: calls[(s,1)] = None
        if (s,2) not in calls: calls[(s,2)] = None

        if calls[(s,1)] is None and calls[(s,2)] is None:
            calls[(s,-1)] = "."
        elif calls[(s,1)] is None:
#            calls[(s,-1)] = str(calls[(s,2)])+"*"
            calls[(s,-1)] = str(calls[(s,2)])
        elif calls[(s,2)] is None:        
#            calls[(s,-1)] = str(calls[(s,1)])+"*"
            calls[(s,-1)] = str(calls[(s,1)])
        elif calls[(s,1)] == calls[(s,2)]:
            calls[(s,-1)] = str(calls[(s,1)])
        else:
            calls[(s,-1)] = "!!"

    #calculate call consistency
    dd2call=set()
    allCall=True
    for s in DD2s:
        if calls[(s,-1)] != ".":
            dd2call.add(calls[(s,-1)])
        else:
            allCall=False
    cons=None
    #all called only if all calls are the same
    if cons is not None:
        allCall = (allCall & cons)
    if len(dd2call)==1: cons=True
    elif len(dd2call)>1: cons=False

#    if rec.is_indel:


    print "\t".join([rec.CHROM,str(rec.POS),str(len(rec.alleles)),rec.var_subtype,str(conc),str(cons),str(allCall)]+
                        [calls[(s,-1)] for s in samples])
    if cons is not None: rec.INFO['DD2CONS']=cons
    if conc is not None: rec.INFO['DD2CONC']=conc
    if allCall is not None: rec.INFO['DD2ALL']=allCall
    if rec.is_snp: rec.INFO['type']='SNP'
    else: rec.INFO['type']='INDEL'
#    rec.INFO['DD2CONS']=str(cons)
#    rec.INFO['DD2CONC']=str(conc)
#    rec.INFO['DD2ALL']=str(allCall)
    vcfout1.write_record(rec)


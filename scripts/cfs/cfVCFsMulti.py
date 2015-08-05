#!/usr/bin/python

import vcf
import sys
import argparse
from string import *
import subprocess
from os.path import basename, splitext
from collections import Counter

parser = argparse.ArgumentParser(description='get STR length differences from disco VCF')

parser.add_argument('-v','--vcf', dest='vcfs', action='append', help='all vcf files for comparison', nargs='+')

parser.add_argument('-o','--out', action="store", dest='outFile', type=str, help='make combined outfile named <outfile> instead of separate files', nargs='?', default=None)
#parser.add_argument('-i','--defaultSize', action="store", dest='defSize', type=int, help='default size for STRs', nargs='?', default=0)


args = parser.parse_args()

vcfs = []
vcfFiles = []
vcfNames = []
vcfSamples = []
sampN = 0


for name in args.vcfs:
    print name
    name = name[0]
    if name.find(':::') > -1:
        vname,fname = name.split(':::')
    else:
        fname = name
        vname = splitext(basename(name))[0]
    vcfFile1 = open(fname,'r')
    reader1=vcf.Reader(vcfFile1)
    vcfNames += [vname]
    vcfs += [reader1]
    vcfFiles += [fname]
    sampN += len(reader1.samples)
    vcfSamples += [reader1.samples] #should be only one

if len(vcfSamples) > len(vcfs):
    print >>sys.stderr, "found more samples than vcfs" 
    print >>sys.stderr, vcfSamples
    print >>sys.stderr, vcfs
    exit(1)  

def _haploidify(alleles):
    if alleles is None: return None
    alls = alleles.split('/')
    
    if len(set(alls)) > 1: print >>sys.stderr,"[warning] could not haploidify "+alleles+", returning "+alls[0] 
    return alls[0]

#merge vcf files:
command = ['vcf-merge']+vcfFiles
print >>sys.stderr, " ".join(command)
tmpvcf = open('merge.'+args.outFile,'w')
ret = subprocess.call(command,stdout=tmpvcf)
tmpvcf.close()
if ret > 0: print >>sys.stderr, "error in "," ".join(command); exit(1)
command = ["bgzip",'merge.'+args.outFile]
ret = subprocess.call(command)
if ret > 0: print >>sys.stderr, "error in "," ".join(command); exit(1)
command = ["tabix",'-pvcf','merge.'+args.outFile+'.gz']
ret = subprocess.call(command)
if ret > 0: print >>sys.stderr, "error in "," ".join(command); exit(1)

print >>sys.stderr,"reading merged VCF: ",'tmp_'+args.outFile+'.gz'

vcfMerge = open('merge.'+args.outFile+'.gz','r')
reader=vcf.Reader(vcfMerge)
rcopy = reader
rcopy.samples = vcfNames 

print >>sys.stderr,"writing results VCF: ",args.outFile
vcfoutF = open(args.outFile,'w')
vcfout=vcf.Writer(vcfoutF,rcopy)

i=0
for rec in reader:
    alleles = [str(a) for a in rec.alleles]
    genos = []
    vartype = "SNP"
    if len(rec.REF) > 1 or len(rec.ALT) > 1:
        vartype="INDEL"
    #print rec.CHROM,rec.POS
    for call in rec.samples:
    #    print call.sample, _haploidify(call.gt_bases)
        allele = _haploidify(call.gt_bases)
        # print alleles, allele, 
        if allele is not None:
        #    print alleles.index(allele)
            genos += [alleles.index(allele)]
        else:
        #    print ""
            genos += [None]
    conc="MISMATCH"    
    if len(alleles) > 2:
        conc="MULTIALLELIC"
    else:
        called = []
        for i in range(0,len(genos)):
            geno = genos[i]
            sample = reader.samples[i]
            if geno==1: #if called as first alternate base
                called += [sample]
#        if (len(called) == len(reader.samples)) and (len(list(set(genos))) == 1):
#            conc="MATCH"+str(len(called))
        conc="MATCH"+str(len(called))
    rec.INFO['CON']=conc
    rec.INFO['MATCHED']=called
    rec.INFO['TYPE']=vartype
    vcfout.write_record(rec)
        
    i+=1
vcfout.close()

def _readVar(reader):
    try:
        rec = reader.next()
    except StopIteration:
        return(None, None, -1, None, None, [], None, True) 
    calls = []
    type = "SNP"
    if len(rec.REF) > 1 or len(rec.ALT) > 1:
        type="INDEL"
    for call in rec.samples:
        calls += [call.gt_type]
    return (rec, rec.CHROM, rec.POS, rec.REF, rec.ALT, calls, type, False)
    

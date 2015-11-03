#!/usr/bin/python

import vcf
import sys
import argparse
from string import *
#from statistics import mean, stdev
from numpy import mean, std
import collections
import tabix
import subprocess
from math import floor, ceil



parser = argparse.ArgumentParser(description='from CHR : POS : 1/0 file, add true/false to VCF')

parser.add_argument('-v','-vcf', action="store", dest='vcfFile', type=str, help='vcf file', nargs='?', default=None)
parser.add_argument('-d','-depthfile', action="store", dest='valueFile', type=str, help='valueFile', nargs='?', default=None)

parser.add_argument('-i','--infoname', action="store", dest='valname', type=str, help='make combined outfile named <outfile> instead of separate files', nargs='?', default="VALUE")
#parser.add_argument('-i','--defaultSize', action="store", dest='defSize', type=int, help='default size for STRs', nargs='?', default=0)


args = parser.parse_args()



vcfFile1 = open(args.vcfFile,'r')
print >>sys.stderr, args.vcfFile, vcfFile1

vcfReader=vcf.Reader(vcfFile1)


_Info = collections.namedtuple('Info', ['id', 'num', 'type', 'desc', 'source', 'version'])
_Filter = collections.namedtuple('Filter', ['id', 'desc'])
_Alt = collections.namedtuple('Alt', ['id', 'desc'])
_Format = collections.namedtuple('Format', ['id', 'num', 'type', 'desc'])

vcfReader.infos['ED'] = _Info(id='ED', num=1, type='Float', desc='Discovar mean edge depth across all samples', source="Discovar", version=1)
vcfReader.infos['EDS'] = _Info(id='EDS', num=1, type='Float', desc='Discovar edge depth std dev across all samples', source="Discovar", version=1)

#vcfReader.infos = infos
for i in vcfReader.infos:
    info = vcfReader.infos[i]
    print "INFO: "+i+" "+str(info)
for format in vcfReader.formats:
    print "FORMAT: "+format
for filt in vcfReader.filters:
    print "FILTER: "+filt




vcfoutF1 = replace(args.vcfFile,'.vcf.gz','.vcf')
vcfoutF1 = replace(vcfoutF1,'.vcf','.REFCALL.vcf')
#print >>sys.stdout, "VCFOUT1",vcfoutF1
vcfoutF1 = open(vcfoutF1,'w')
vcfout1=vcf.Writer(vcfoutF1,vcfReader)

valfile = tabix.open(args.valueFile)

vals = dict()
dns = []

outGenoCalldata = vcf.model.make_calldata_tuple("GT")

#header = subprocess.check_output(["zcat",args.valueFile,"|","head","-n","1"],shell=True)
#dns = header[2:]

zcat = subprocess.Popen(['zcat',args.valueFile], stdout=subprocess.PIPE)
header = subprocess.check_output(['head','-n','1'], stdin=zcat.stdout)
print header
dns=split(header)[2:]
print dns

def _getDepthForChrom(chrom):
    vals = _getDepthForRegion(chrom,0,10000000)
    return vals

def _getDepthForRegion(chrom,st,en):
    print >>sys.stderr, "parsing depth file for "+chrom+":"+str(st)+"-"+str(en)
    vals = dict()

    for rs in valfile.query(chrom,st,en):
#        rs = vline.split()
#        print >>sys.stderr, rs

        c,p = rs[:2]
        ds = rs[2:]
        
#        if c == "CHR":
#            dns = ds
#        else:
        ds = map(int,ds)
        for i in range(0,len(ds)):
            dn = dns[i]
            vals[(c,int(p),dn)] = ds[i]
    #print vals.keys()[:10]
    return vals

thisChrom=""
thisMax=0
blockSize=1e5

print >>sys.stderr, "parsing VCF file"
for rec in vcfReader:
    #if rec.is_indel:

    if thisChrom != rec.CHROM:
        newMax = int(ceil(rec.POS/blockSize)*blockSize)
        vals = _getDepthForRegion(rec.CHROM,0,newMax)
        thisChrom = rec.CHROM
        thisMax = newMax
    elif rec.POS > thisMax:
        newMax = int(ceil(rec.POS/blockSize)*blockSize)
        vals = _getDepthForRegion(rec.CHROM,thisMax,newMax)
        thisChrom = rec.CHROM
        thisMax = newMax

    samples = rec.samples

    


    allDepths = []
    for i in range(0,len(samples)):
        call = samples[i]
        depth = -1
        if (rec.CHROM,rec.POS,call.sample) in vals:
            depth = vals[(rec.CHROM,rec.POS,call.sample)]
            allDepths += [depth]
            
#            rec.INFO[args.valname]=vals[(rec.CHROM,rec.POS)]
        if call.gt_type is None and depth > 0:
#        if call.gt_type is None and (rec.CHROM,rec.POS,call.sample) in vals:
            #make new tuple with ref genotype
            cd = outGenoCalldata("0/0")
            #make array of required vals (record, samplename, calldata)
            # samples += [vcf.model._Call(rec,call.sample,cd)] 
            samples[i] = vcf.model._Call(rec,call.sample,cd) 
#            print >>sys.stderr, '.',
        elif call.gt_type is None and depth == -1:
            # print >>sys.stderr, "warning: ",rec.CHROM+":"+str(rec.POS)+" "+call.sample+" not found in valfile and GT is null"
#            print >>sys.stderr, '+',
            pass
            #leave as null
            # cd = outGenoCalldata("./.")
            # samples[i] = [vcf.model._Call(rec,call.sample,cd)] 
        
    rec.samples = samples
    rec.INFO['ED'] = mean(allDepths)
    rec.INFO['EDS']= round(std(allDepths),3)
    vcfout1.write_record(rec)


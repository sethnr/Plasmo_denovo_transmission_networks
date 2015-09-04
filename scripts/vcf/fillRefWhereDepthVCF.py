#!/usr/bin/python

import vcf
import sys
import argparse
from string import *
#from statistics import mean, stdev
from numpy import mean, std
import collections 
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

valfile = open(args.valueFile,'r')

vals = dict()
dns = []

outGenoCalldata = vcf.model.make_calldata_tuple("GT")

print >>sys.stderr, "parsing depth file"
for vline in valfile:
    rs = vline.split()
    c,p = rs[:2]
    ds = rs[2:]

#    print c,"!",p
    if c == "CHR":
        dns = ds
    else:
        ds = map(int,ds)
    #if p == 451480: print c,p,v,bool(v) 
        for i in range(0,len(ds)):
            dn = dns[i]
            vals[(c,int(p),dn)] = ds[i]


print >>sys.stderr, "parsing VCF file"
for rec in vcfReader:
    #if rec.is_indel:

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
        elif call.gt_type is None and depth == -1:
           print >>sys.stderr, "warning: ",rec.CHROM+":"+str(rec.POS)+" "+call.sample+" not found in valfile and GT is null"
            #leave as null
            # cd = outGenoCalldata("./.")
            # samples[i] = [vcf.model._Call(rec,call.sample,cd)] 
        
    rec.samples = samples
    rec.INFO['ED'] = mean(allDepths)
    rec.INFO['EDS']= round(std(allDepths),3)
    vcfout1.write_record(rec)


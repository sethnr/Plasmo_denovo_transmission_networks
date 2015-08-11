#!/usr/bin/python

import vcf
import sys
import argparse
from string import *
import collections
parser = argparse.ArgumentParser(description='get STR length differences from disco VCF')

parser.add_argument('-v','--vcf', action="store", dest='vcfFile1', type=str, help='vcfFile1', nargs='?', default=None)
parser.add_argument('-s','--samples', action="append", dest='samples', type=str, help='samples to calc across', nargs='?', default=None)
parser.add_argument('-p','--prefixVars', action="store", dest='prefix', type=str, help='samples to calc across', nargs='?', default="")

parser.add_argument('-o','--out', action="store", dest='outFile', type=str, help='make combined outfile named <outfile> instead of separate files', nargs='?', default=None)
#parser.add_argument('-i','--defaultSize', action="store", dest='defSize', type=int, help='default size for STRs', nargs='?', default=0)


args = parser.parse_args()

vcfFile1 = open(args.vcfFile1,'r')
reader1=vcf.Reader(vcfFile1)


_Info = collections.namedtuple('Info', ['id', 'num', 'type', 'desc', 'source', 'version'])
_Filter = collections.namedtuple('Filter', ['id', 'desc'])
_Alt = collections.namedtuple('Alt', ['id', 'desc'])
_Format = collections.namedtuple('Format', ['id', 'num', 'type', 'desc'])

reader1.infos[args.prefix+'MISS'] = _Info(id=args.prefix+'MISS', num=1, type='Float', desc=args.prefix+' missingness', source="checkLaneConcordance.py", version=1)
reader1.filters[args.prefix+'LaneErr'] = _Filter(id=args.prefix+'LaneErr', desc=args.prefix+' lanes not concordant')


reader1.infos['ED'] = _Info(id='ED', num=1, type='Float', desc='Discovar mean edge depth across all samples', source="Discovar", version=1)
reader1.infos['EDS'] = _Info(id='EDS', num=1, type='Float', desc='Discovar edge depth std dev across all samples', source="Discovar", version=1)


# for i in reader1.infos:
#    info = reader1.infos[i]
#    print "INFO: "+i+" "+str(info)
#for format in reader1.formats:
#    print "FORMAT: "+format
#for filt in reader1.filters:
#    print "FILTER: "+filt


vcfoutF1 = replace(args.vcfFile1,'.vcf','.'+args.prefix+'LCHK.vcf')
vcfoutF1 = replace(vcfoutF1,'.vcf.gz','.'+args.prefix+'LCHK.vcf')
print >>sys.stdout, "VCFOUT1",vcfoutF1
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

print "\t".join(["chr","pos","alleles","type"]+samples)
for rec in reader1:
    calls = dict()
    for call in rec.samples:
        sample, lane = call.sample.split('_')
        calls[(sample,int(lane))] = call.gt_type
#    samples = list(set([s for s,l in calls]))
    laneError = 0
    missing = 0.0
    for s in samples:
        if (s,1) not in calls: calls[(s,1)] = None
        if (s,2) not in calls: calls[(s,2)] = None

        if calls[(s,1)] is None and calls[(s,2)] is None:
            calls[(s,-1)] = "."
            missing += 2
        elif calls[(s,1)] is None:
            calls[(s,-1)] = str(calls[(s,2)])+"*"
            missing += 1
        elif calls[(s,2)] is None:        
            calls[(s,-1)] = str(calls[(s,1)])+"*"
            missing += 1
        elif calls[(s,1)] == calls[(s,2)]:
            calls[(s,-1)] = str(calls[(s,1)])
        else:
            calls[(s,-1)] = "!!"
            laneError +=1
#    if rec.is_indel:
#        print "\t".join([rec.CHROM,str(rec.POS),str(len(rec.alleles)),rec.var_subtype]+
#                        [calls[(s,-1)] for s in samples])

    rec.INFO[args.prefix+'MISS']=round(missing/len(rec.samples),4)
    if laneError > 0:
#        print len(rec.FILTER)
        if len(rec.FILTER) == -1:
            rec.FILTER[1]=args.prefix+'LaneErr'
        else:
            rec.FILTER += [args.prefix+'LaneErr']

    if len(rec.ALT) > 50:
        altp = reader1.formats.keys().index('ALTP')
        print "\t".join([rec.CHROM,str(rec.POS),str(len(rec.ALT)),",".join(rec.FILTER)])
#        for call in rec.samples:
#            print call.data[altp]
    else:
        vcfout1.write_record(rec)

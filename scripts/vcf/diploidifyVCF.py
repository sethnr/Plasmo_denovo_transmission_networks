#!/usr/bin/python

import vcf
import sys
import argparse
from string import *
import collections
import copy

parser = argparse.ArgumentParser(description='get STR length differences from disco VCF')

parser.add_argument('-v','--vcf', action="store", dest='vcfFile1', type=str, help='vcfFile1', nargs='?', default=None)

args = parser.parse_args()

vcfFile1 = open(args.vcfFile1,'r')
reader1=vcf.Reader(vcfFile1)

vcfoutF1 = replace(args.vcfFile1,'.vcf','.DIPLOID.vcf')
vcfoutF1 = replace(vcfoutF1,'.vcf.gz','.vcf')
print >>sys.stdout, "VCFOUT1",vcfoutF1
vcfoutF1 = open(vcfoutF1,'w')

reader2 = copy.copy(reader1)

vcfout1=vcf.Writer(vcfoutF1,reader2)

makeCallData = vcf.model.make_calldata_tuple(("GT","ALTP","REFP","GP"))

print >>sys.stderr, reader1.formats


GTindex=None
for rec in reader1:
#    calls = dict()
#    rec.samples = None
    newSamples = list()
    GTindex=rec.FORMAT.split(":").index("GT")
    format = rec.FORMAT.split(":")
    makeCallData = vcf.model.make_calldata_tuple(tuple(rec.FORMAT.split(":")))
    
    for call in rec.samples:
#        print call.data
#        print rec.FORMAT
        callData= list(call.data)
        dipGT = callData[GTindex]
        
        if dipGT is not None:
            dipAlleles = dipGT.split('/')
            if(len(dipAlleles)==1):
                callData[GTindex]=str(dipAlleles[0])+"/"+str(dipAlleles[0])
        
        newCallString = ""

        for i in range(0,len(callData)):
            if callData[i] is None: callData[i]='.'
            if type(callData[i]) is list: callData[i]=",".join(map(str,callData[i]))
            newCallString += format[i]+'="'+str(callData[i])+'", '
        newCallData = eval("makeCallData("+newCallString+")")
        print call.data
        print newCallData
        newCall = vcf.model._Call(rec,call.sample,newCallData) 
        newSamples += [newCall]
    print rec.samples
    rec.samples=newSamples
    print rec.samples
    
#makeCallData(namedtuple(callData))
        

        #print >> sys.stderr,  call.sample
#        if sampleKey[call.sample] not in newSampleNames:
#            newSampleNames += [sampleKey[call.sample]]
#            call.sample=sampleKey[call.sample]
#            newSamples += [call]
    #print >>sys.stderr, len(rec.samples), len(newSamples)
#    rec.samples = newSamples
    vcfout1.write_record(rec)

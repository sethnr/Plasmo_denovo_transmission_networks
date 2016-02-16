#!/usr/bin/python

import vcf
import sys
import argparse
from string import *
import collections
import copy

parser = argparse.ArgumentParser(description='mark singletons in file')

parser.add_argument('-v','--vcf', action="store", dest='vcfFile1', type=str, help='vcfFile1', nargs='?', default=None)
parser.add_argument('-R','--keepRefCalls', action="store_true", dest='keepRefCalls', help='do not remove calls in which only ref is called', default=False)

args = parser.parse_args()

vcfFile1 = open(args.vcfFile1,'r')
reader1=vcf.Reader(vcfFile1)

vcfoutF1 = replace(args.vcfFile1,'.vcf','.MKSNGL.vcf')
vcfoutF1 = replace(vcfoutF1,'.vcf.gz','.vcf')
print >>sys.stdout, "VCFOUT1",vcfoutF1
vcfoutF1 = open(vcfoutF1,'w')


_Filter = collections.namedtuple('Filter', ['id', 'desc'])
reader1.filters['Singleton'] = _Filter(id='Singleton', desc='only one minor variant at locus')
reader1.filters['Doubleton'] = _Filter(id='Doubleton', desc='only one minor variant at locus')
reader1.filters['Tripleton'] = _Filter(id='Tripleton', desc='only one minor variant at locus')

reader2 = copy.copy(reader1)

vcfout1=vcf.Writer(vcfoutF1,reader2)

#makeCallData = vcf.model.make_calldata_tuple(("GT","ALTP","REFP","GP"))

print >>sys.stderr, reader1.formats

doubletons = dict()
for x in reader1.samples:
    doubletons[x] = dict()
    for y in reader1.samples:
        doubletons[x][y]=0

GTindex=None
for rec in reader1:
#    calls = dict()
#    rec.samples = None
    newSamples = list()
    GTindex=rec.FORMAT.split(":").index("GT")
    GPindex=rec.FORMAT.split(":").index("GP")
    format = rec.FORMAT.split(":")
    makeCallData = vcf.model.make_calldata_tuple(tuple(rec.FORMAT.split(":")))
    
    callCount = dict()
#    print rec.CHROM, rec.POS, len(rec.ALT),rec.is_indel,
    for call in rec.samples:
        callData= list(call.data)
        GT = callData[GTindex]
        if GT is not None:
#            print "\t"+GT,
            if GT not in callCount:
                callCount[GT]=0
            callCount[GT] +=1
        else:
            pass
#            print "\t.",

    if len(callCount)==1:
        #REF CALL ONLY
#        print "\tREFCALL - "+str(callCount.keys()[0])
        if callCount.keys()[0] == 0:
            rec.FILTER += ["RefCallOnly"]
        if args.keepRefCalls:
            vcfout1.write_record(rec)
        
    else:
        sortArr = [(callCount[c], c) for c in callCount] 
        sortArr.sort(reverse=True)
        sortKeys = [c for (n,c) in sortArr]
        nonRef = 0
        for call in sortKeys[1:]:
            if callCount[call]==1:
#                if nonRef==1:
#                    rec.FILTER += ["Singleton"]
#                print "\tSINGLETON"
                pass
            if callCount[call]==2:
                rec.FILTER += ["Doubleton"]
#                print "\tDOUBLETON"
            
                shareDoub = list()
                for geno in rec.samples:
                    callData= list(geno.data)
                    GT = callData[GTindex]
                    if GT == call:
                        shareDoub += [geno.sample]
                doubletons[shareDoub[0]][shareDoub[1]] += 1
                doubletons[shareDoub[1]][shareDoub[0]] += 1
#        if nonRef==3:
#            rec.FILTER += ["Tripleton"]
#            print "\tTRIPLETON"
            else:
                pass
#                print "\t"+str(nonRef)
#        vcfout1.write_record(rec)
    
#makeCallData(namedtuple(callData))

for y in reader1.samples:
    print >>sys.stdout,"\t",y,
print >>sys.stdout,""
for x in reader1.samples:
    print >>sys.stdout,x,
    for y in reader1.samples:
        print >>sys.stdout,"\t",doubletons[x][y],
    print >>sys.stdout,""
        

        #print >> sys.stderr,  call.sample
#        if sampleKey[call.sample] not in newSampleNames:
#            newSampleNames += [sampleKey[call.sample]]
#            call.sample=sampleKey[call.sample]
#            newSamples += [call]
    #print >>sys.stderr, len(rec.samples), len(newSamples)
#    rec.samples = newSamples

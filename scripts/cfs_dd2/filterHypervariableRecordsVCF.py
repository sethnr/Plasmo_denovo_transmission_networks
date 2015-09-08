#!/usr/bin/python

import vcf
import sys
import argparse
from string import *
import collections
import numpy as np

parser = argparse.ArgumentParser(description='get STR length differences from disco VCF')

parser.add_argument('-v','--vcf', action="store", dest='vcfFile1', type=str, help='vcfFile1', nargs='?', default=None)
parser.add_argument('-p','--minProb', action="store", dest='minProb', type=float, help='minimum probability above which we retain allele (p=0 -> retain all positive)', nargs='?', default=0.0)
parser.add_argument('-F','--noFilter', action="store_false", dest='filter', help='do not add filter if MAF > 1', default=True)


args = parser.parse_args()

vcfFile1 = open(args.vcfFile1,'r')
reader1=vcf.Reader(vcfFile1)

_Info = collections.namedtuple('Info', ['id', 'num', 'type', 'desc', 'source', 'version'])
_Filter = collections.namedtuple('Filter', ['id', 'desc'])
_Alt = collections.namedtuple('Alt', ['id', 'desc'])
_Format = collections.namedtuple('Format', ['id', 'num', 'type', 'desc'])

reader1.infos['ANO'] = _Info(id='ANO', num=1, type='Integer', desc='potential alleles predicted by discovar (prior to filtering)', source="filterHypervariables", version=1)
if args.filter:
    reader1.filters['Hypervariable'] = _Filter(id='Hypervariable', desc='more than 50 alleles in vcf')


vcfoutF1 = replace(args.vcfFile1,'.vcf','.HYPF.vcf')
vcfoutF1 = replace(vcfoutF1,'.vcf.gz','.vcf')
vcfoutF1 = open(vcfoutF1,'w')
vcfout1=vcf.Writer(vcfoutF1,reader1)

samples = reader1.samples

#maxAlleles=len(samples)
maxAlleles=50


makeCallData = vcf.model.make_calldata_tuple(("GT","ALTP","REFP","GP"))


for rec in reader1:
    calls = dict()
    altCount = len(rec.ALT)
    
    if len(rec.ALT) > maxAlleles:
        #get index of altp and refp sets
        format = rec.FORMAT.split(":")
        altpi = format.index('ALTP')
        refpi = format.index('REFP')
        gti = format.index('GT')
        gpi = format.index('GP')
#        print altp
#        print refp
#        print "\t".join([rec.CHROM,str(rec.POS),str(len(rec.ALT)),",".join(rec.FILTER)])

        #MAKE MATRIX OF ALLELE PROBABILITIES
        aProbs = np.zeros((altCount,len(samples)),dtype=np.float)
        si=0
        
        genotyped = [False]*altCount
        for call in rec.samples:
            altprobs = call.data[altpi]
            refprob = call.data[refpi]
            # set genotyped = true if has a genotype (ignore phred probs)
            
            gt = call.data[gti]
            if gt is not None:
                a1,a2 = gt.split('/')
                a1 = int(a1)
                a2 = int(a2)
                if a1 > 0:
#                    print a1
                    genotyped[a1-1]=True
                if a2 > 0:
                    genotyped[a2-1]=True
            
            if refprob is None:
                probs = [0]*altCount
            else:
                probs = altprobs
            #replace None vals with zero
            probs = [p if p is not None else 0 for p in probs]
            probs = map(float,probs)
            aProbs[:,si] = probs
            si +=1
        overProb = np.greater(np.sum(aProbs,1),args.minProb)

#        print rec.POS
#        print "PROBS: " ,
#        print np.flatnonzero(overProb)
#        print "GENOS: ",
#        print np.array(genotyped)
        passIndex = np.flatnonzero(np.logical_or(overProb,genotyped))
#        print passIndex
        allLookup = dict()
        allLookup[0]=0  #ref doesn't change position

        newALT = []
        ii = 1
        for i in passIndex.tolist():
            i = i+1
            #rec.ALT[i]
            allLookup[i]=ii
            ii+=1
#            print str(i)+"-"+str(ii)
            newALT += [rec.ALT[i-1]]
        if newALT == []:
            newALT += [rec.ALT[0]]
            
        newCalls = []
        for call in rec.samples:
            #REPLACE WITH NEW ALLELE VALUES
            #NEED TO MAKE NEW CALLDATA OBJECT?!
            #print "GT bases: ", call.gt_bases 
            if call.gt_bases is None:
                newgt = './.'
            else:
                oldgt = call.data[gti]
                a1,a2 = oldgt.split('/')
                a1 = int(a1)
                a2 = int(a2)
            #    print >>sys.stderr, "GT:", oldgt
                try:
                    newgt = str(allLookup[a1])+'/'+str(allLookup[a2])
                except KeyError:
                    if a1 in allLookup:
                        newgt = str(allLookup[a1])+'/'+str(allLookup[a1])
                    elif a2 in allLookup:
                        newgt = str(allLookup[a2])+'/'+str(allLookup[a2])
                    else:
                        newgt = './.'
                    print >>sys.stderr, "KeyError: "+str(a1)+"/"+str(a2)+" - "+rec.CHROM+":"+str(rec.POS)+" allele prob not above "+str(args.minProb)+" making call "+newgt
#            print call.data

            altP = call.data[altpi]
            refP = call.data[refpi]
            gp = call.data[gpi]
            
#            print len(rec.ALT)
#            print len(newALT)
            if altP is not None:
#                print altP, altpi
#                print "ALTP LEN: "+str(len(altP))
                newAltP = []
                for i in allLookup.keys():
#                    print i, allLookup[i], altP[allLookup[i]],  altP[allLookup[i]-1]
                    newAltP += [altP[allLookup[i]]]
                altP = newAltP
            if gp is not None:
#                print gp, gpi
#                print "GP LEN: "+str(len(gp))
                gp = None   #remove to make file more readable - frankly have no idea what the numbers refer to anyway
            newCallData = makeCallData(GT=newgt,ALTP=altP,REFP=refP,GP=gp)
            newCall = vcf.model._Call(rec,call.sample,newCallData) 
            newCalls += [newCall]
        rec.ALT=newALT
        rec.samples = newCalls
                               
#        print np.array(rec.ALT)[np.where(overProb[1:])]
        if args.filter:
            if len(rec.FILTER) == -1:
                rec.FILTER[1]='Hypervariable'
            else:
                rec.FILTER += ['Hypervariable']
        
        
    
    rec.INFO['ANO'] = altCount
    vcfout1.write_record(rec)


#!/usr/bin/python

import sys
import argparse
from string import *
from subprocess import call, check_output
import os.path as path
from math import floor,ceil


import numpy as np
import vcf
#from pysam import TabixFile
import tabix

parser = argparse.ArgumentParser(description='get STR length differences from disco VCF')

parser.add_argument('-v','--vcf', action="store", dest='vcfFile', type=str, help='vcf file', nargs='?', default=None)
parser.add_argument('-d','--depth', action="store", dest='depth', type=str, help='output of samtools depth', nargs='?', default=None)

parser.add_argument('-o','--out', action="store", dest='outFile', type=str, help='outFile', nargs='?', default=None)
parser.add_argument('-b','--blockSize', action="store", dest='blocksize', type=int, help='default blocksize for summary (kb)', nargs='?', default=1)
parser.add_argument('-c','--minCov', action="store", dest='minCovDepth', type=int, help='minimum coverage to consider something covered', nargs='?', default=5)

parser.add_argument('-r','--region', action="store", dest='regions', type=str, help='regions over which to summarize', nargs='+', default=None)

args = parser.parse_args()
print >>sys.stderr, args.vcfFile, args.depth


outfile = sys.stdout
if args.outFile is not None:
    outfile = open(args.outFile,'w')

def _getNoSNPs(genos,blockSize):
    return genos.shape[2] / blockSize


def _getHaploidPi(matrix):
    pass


if not path.isfile(args.vcfFile+".tbi"):
    #make tabix index for vcf if not tabixed
    if path.splitext(args.vcfFile) != '.gz':
        call(["bgzip",args.vcfFile])
        args.vcfFile += '.gz'
    call(["tabix","-p","vcf",args.vcfFile])
vcfFile = open(args.vcfFile,'r')
reader=vcf.Reader(vcfFile)


if not path.isfile(args.depth+".tbi"):
    #make tabix index for depth if not tabixed
    if path.splitext(args.depth) != '.gz':
        call(["bgzip",args.depth])
        args.depth += '.gz'
    call(["tabix","-b2","-e2",args.depth])
#depth = sam.TabixFile(args.depth)
depth = tabix.open(args.depth)

if args.regions is None:
    #get regions from vcffile
    pass

##POSSIBLE VALUES FOR CONCORDANCE AND MATCHED IN VCF
contypes = ['MATCH0','MATCH1','MATCH2','MATCH3','MULTIALLELIC']
matchtypes = ['NUCMER','DISCO','HAPLO',
              'NUCMER,DISCO','NUCMER,HAPLO','DISCO,HAPLO',
              'NUCMER,DISCO,HAPLO','']

#blocksize in kb => block in bases
block=args.blocksize * 1000
for region in args.regions:
    (chrom, locus) = region.split(":")
    (regst, regen) = locus.split("-")
    regst = int(floor(int(regst) / block)*block)
    regen = int(ceil(int(regen) / block)*block)
    for st in range(regst,regen,block):
        en = st+block

        ##################
        # PROCESS DEPTH VALUES
        #################
        #get depth vals from (tabixed) depth file
        depths = None
        i=-1
        for bdepths in depth.query(chrom, st, en):
            i+=1
##            bdepths = row.split()
            bdepths = bdepths[2:]
#            print bdepths, len(bdepths),block
            #initialise empty zero arrays
            if depths is None:
                #initialise depths
                depths = np.zeros((block,len(bdepths)))
                if st == regst:
                    #print header:
                    print >>outfile,"\t".join(["chrom","start","end"]+
                                              ["cov"+str(d+1) for d in range(0,len(bdepths))]+
                                              ["d"+str(d+1) for d in range(0,len(bdepths))]+
                                              ["dcov"+str(d+1) for d in range(0,len(bdepths))]+
                                              ["SNPS"]+["S-"+c for c in contypes] +
                                              ["S-"+m for m in matchtypes] +
                                              ["INDELS"]+["I-"+c for c in contypes] +
                                              ["I-"+m for m in matchtypes]
                                              )
            for j in range(0,len(bdepths)):
#                print >>sys.stderr,"I:",i," J:",j
                bdepth = float(bdepths[j])
                depths[(i,j)] = bdepth
##        print >>sys.stderr, depths
        mean = np.mean(depths,axis=0)

#        depths[depths < args.minCovDepth] = np.nan
        #print >>sys.stderr, st,en, depths, args.minCovDepth, sum(depths,0)
#        if np.all(np.isnan(depths)):
#            noCovered = 0
#        else:
            #noCovered = np.sum(depths >= args.minCovDepth,0)
#            noCovered = np.sum(depths[np.isfinite(depths)] >= args.minCovDepth,0)
        noCovered = np.sum(depths >= args.minCovDepth,0)
        #print >>sys.stderr, noCovered
        depths[depths < args.minCovDepth] = np.nan
        meanCovered = np.nanmean(depths,axis=0)
        meanCovered[np.isnan(meanCovered)] = 0
        
        print >>outfile, "\t".join([chrom, str(st), str(en)] +
                                   [str(i) for i in noCovered.tolist()+mean.tolist()+meanCovered.tolist()]),

        ##################
        # PROCESS GENOTYPES
        #################
        blockreader = reader.fetch(chrom, st, en)
        vars = 0
        #setup null count matrix

        typecon = dict()
        typematch = dict()
        total=dict()
        total['SNP']=0
        total['INDEL']=0

        for vartype in ['SNP','INDEL']:
#            for con in ['MATCH','MISMATCH','PRIVATE']:
            for con in contypes:
                typecon[(vartype,con)] = 0
            for match in matchtypes:
                typematch[(vartype,match)] = 0

        for rec in blockreader:
            vars += 1
            vartype = rec.INFO['TYPE']
            total[vartype] += 1
            con = rec.INFO['CON'][0]
            match = ','.join(rec.INFO['MATCHED'])
        
#            if (vartype, con) not in typecon:
#                typecon[(vartype,con)] = 0
            typecon[(vartype,con)] +=1
            typematch[(vartype,match)] +=1
            
        for vartype in ['SNP','INDEL']:
#            total = 0
#            for con in ['MATCH','MISMATCH','PRIVATE']:
#            for con in contypes:
#                total += typecon[vartype,con]
            print >>outfile, "\t"+str(total[vartype]),
#            for con in ['MATCH','MISMATCH','PRIVATE']:
            for con in contypes:
                if total[vartype] > 0:
                    pc = typecon[vartype,con]/float(total[vartype])
                else:
                    pc=float(0)
                print >>outfile,"\t"+str(pc),
            for match in matchtypes:
                if total[vartype] > 0:
                    pc = typematch[vartype,match]/float(total[vartype])
                else:
                    pc=float(0)
                print >>outfile,"\t"+str(pc),
        print ""

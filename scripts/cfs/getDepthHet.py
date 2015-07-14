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

    
vcfFile = open(args.vcfFile,'r')
reader=vcf.Reader(vcfFile)

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
                    print >>outfile,"\t".join(["#chrom","start","end"]+
                                          ["d"+str(d+1) for d in range(0,len(bdepths))]+
                                          ["SNPS","Smatch","Smiss","Spriv",
                                          "INDELS","Imatch","Imiss","Ipriv"])
            for j in range(0,len(bdepths)):
#                print i,j
                bdepth = float(bdepths[j])
                depths[(i,j)] = bdepth
##        print >>sys.stderr, depths
        mean = np.mean(depths,axis=0)

        depths[depths < args.minCovDepth] = np.nan
        noCovered = np.sum(depths >= args.minCovDepth,0)
        meanCovered = np.nanmean(depths,axis=0)
        meanCovered[np.isnan(meanCovered)] = 0
        
        print >>outfile, chrom, st, en, "\t".join([str(i) for i in noCovered.tolist()+mean.tolist()+meanCovered.tolist()]),

        ##################
        # PROCESS GENOTYPES
        #################
        blockreader = reader.fetch(chrom, st, en)
        vars = 0
        #setup null count matrix
        typecon = dict()
        for type in ['SNP','INDEL']:
            for con in ['MATCH','MISMATCH','PRIVATE']:
                typecon[(type,con)] = 0

        for rec in blockreader:
            vars += 1
            type = rec.INFO['TYPE'][0]
            con = rec.INFO['CON'][0]
            if (type, con) not in typecon:
                typecon[(type,con)] = 0
            typecon[(type,con)] +=1
        
        for type in ['SNP','INDEL']:
            total = 0
            for con in ['MATCH','MISMATCH','PRIVATE']:
                total += typecon[type,con]
            print >>outfile, total,
            for con in ['MATCH','MISMATCH','PRIVATE']:
                if total > 0:
                    pc = typecon[type,con]/float(total)
                else:
                    pc=float(0)
                print >>outfile,pc,
        print ""

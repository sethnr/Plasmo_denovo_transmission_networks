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

        #GET DEPTH VALUES
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
                depths = np.zeros((block,len(bdepths)))

            for j in range(0,len(bdepths)):
                bdepth = int(bdepths[j])
                depths[(i,j)] = bdepth
##        print >>sys.stderr, depths
        


        #GET GENOTYPES
        blockreader = reader.fetch(chrom, st, en)
        for rec in blockreader:
            pass
    

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
parser.add_argument('-i','--infoname', action="store", dest='valname', type=str, help='make outfile named <vcfTo>.<infoname>.vcf. ', nargs='+', default=None)

parser.add_argument('-b','--blockSize', action="store", dest='blocksize', type=int, help='default blocksize for summary (kb)', nargs='?', default=1)

parser.add_argument('-r','--region', action="store", dest='regions', type=str, help='regions over which to summarize', nargs='+', default=None)

args = parser.parse_args()
outfile = sys.stdout

if not path.isfile(args.vcfFile+".tbi"):
    #make tabix index for vcf if not tabixed
    if path.splitext(args.vcfFile) != '.gz':
        call(["bgzip",args.vcfFile])
        args.vcfFile += '.gz'
    call(["tabix","-p","vcf",args.vcfFile])
vcfFile = open(args.vcfFile,'r')
reader=vcf.Reader(vcfFile)


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
            for con in contypes:
                typecon[(vartype,con)] = 0
            for match in matchtypes:
                typematch[(vartype,match)] = 0

        for rec in blockreader:
            vars += 1
            match = ','.join(rec.INFO['MATCHED'])
            
            typecon[(vartype,con)] +=1
            typematch[(vartype,match)] +=1
            
        for vartype in ['SNP','INDEL']:
            print >>outfile, "\t"+str(total[vartype]),
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

#!/usr/bin/python

import sys
import argparse
from string import *
from subprocess import call, check_output
import os.path as path
from math import floor,ceil
import bz2

import numpy as np
import vcf
#from pysam import TabixFile
import tabix

parser = argparse.ArgumentParser(description='get STR length differences from disco VCF')

parser.add_argument('-v','--vcf', action="store", dest='vcfFile', type=str, help='vcf file', nargs='?', default=None)
parser.add_argument('-o','--outname', action="store", dest='valname', type=str, help='make outfile named <vcfTo>.<infoname>.vcf. ', nargs='+', default=None)

#parser.add_argument('-b','--blockSize', action="store", dest='blocksize', type=int, help='default blocksize for summary (kb)', nargs='?', default=1)

parser.add_argument('-r','--region', action="store", dest='regions', type=str, help='regions over which to summarize', nargs='+', default=None)

parser.add_argument('-i','--infofield', action="append", dest='info', type=str, help='invo values which should be summarized', default=None)

args = parser.parse_args()
outfile = sys.stdout

## if not path.isfile(args.vcfFile+".tbi"):
##    #make tabix index for vcf if not tabixed
##    if path.splitext(args.vcfFile) != '.gz':
##        call(["bgzip",'-f',args.vcfFile])
##        args.vcfFile += '.gz'
##    call(["tabix","-p","vcf",args.vcfFile])

#SCAN THROUGH TO GET ALL FILTERS
vcfFile = open(args.vcfFile,'r')
reader=vcf.Reader(vcfFile)
filters = list()
for rec in reader:
    filters += rec.FILTER
filters = list(set(filters))

#print "ALL FILTERS"
#print filters
vcfFile.close()


#REOPEN SAME FILE
vcfFile = open(args.vcfFile,'r')
reader=vcf.Reader(vcfFile)

def stringComplex(mystring):
    comp = bz2.compress(mystring)
    bcomplex = float(len(mystring))/len(comp)
    return bcomplex

def _basecomp(allele):
    Apc = sum([1 for i in allele if i in ('a','A')])/float(len(allele))
    Tpc = sum([1 for i in allele if i in ('t','T')])/float(len(allele))
    Cpc = sum([1 for i in allele if i in ('c','C')])/float(len(allele))
    Gpc = sum([1 for i in allele if i in ('g','G')])/float(len(allele))
    basePC = ":".join([str(round(i,2)) for i in [Apc,Tpc,Cpc,Gpc]])
    return (Apc,Tpc,Cpc,Gpc)

#print header:
print >>sys.stdout, "\t".join(
            ["chr",
            "pos",
            "vartype",
            #                    '.'.join(alleles),
            "varlen",
            "vcomplex",
            "STRtype",
            "INDtype",
            "DUST",
            "STRP",
            "STRE",
            "alleles",
            "call_rate",
            "nuc_div"]+
            filters+
            ["coding",
            "consequence",
            'DD2FDKMAF'])
 


for rec in reader:
    if rec.is_indel:
        vartype='INDEL'
    elif rec.is_snp:
        vartype='SNP'
    else:
        print >>sys.stderr, rec.CHROM+":"+str(rec.POS)+"-"+rec.REF+'/'.join(rec.ALT)+" is neither SNP or INDEL"
    STRtype=""
    STRP=""
    STRE=""
    if 'STR' in rec.INFO:
        baseComp = rec.INFO['STRatcg'][0]
#        print baseComp
        fa, fc, ft, fg = map(float,baseComp.split(':'))
        STRtype="STR"
        if fa == 1:
            STRtype="polyA"
        elif ft == 1:
            STRtype="polyA"
        elif fa >= 0.4 and ft >= 0.4:
            STRtype="TArep"
        elif fc == 0 and fg==0:
            STRtype="TArich"
        
        STRP=rec.INFO['STRP'][0]
        STRE=rec.INFO['STRE'][0]
    dust=0
    if 'DUST' in rec.INFO:
        dust=1

    consequence=""
    coding="intergenic"
    if 'ANN' in rec.INFO:        
        coding="coding"
        snpEff = rec.INFO['ANN'][0]
#        print snpEff
        Allele, Annotation, Annotation_Impact, Gene_Name, Gene_ID, Feature_Type, Feature_ID, Transcript_BioType, Rank, HGVSc, HGVSp, cDNAposlen, CDSposlen, AAposlen, Distance, Errors = snpEff.split('|')
        consequence=Annotation
    
        
    INDtype = ""
    #get indel complexity
    varlen = 1
    bcomplex=0
    alleles = [str(i) for i in rec.alleles]
    if rec.is_indel:        
        reflen = len(alleles[0])
        maxallele = max(alleles[1:], key=len)
        altlen = len(maxallele)
        comp = bz2.compress(maxallele)
        bcomplex = round(stringComplex(maxallele),3)
        varlen=altlen-reflen

        if varlen < 0:
            delseq = alleles[0]
            basecomp = _basecomp(delseq[1:])
        elif varlen >0:
            inseq = maxallele
            basecomp = _basecomp(inseq[1:])

        fa,ft,fc,fg = basecomp
        if fa == 1:
            INDtype="polyA"
        elif ft == 1:
            INDtype="polyA"
        elif fa >= 0.4 and ft >= 0.4:
            INDtype="TArep"
        elif fc == 0 and fg==0:
            INDtype="TArich"
  
            
    #else:
    #    basecomp=_basecomp(alleles[1])

    

    #get no of alleles
    alleleCount=len(alleles)
    if 'ANO' in rec.INFO:
        alleleCount=int(rec.INFO['ANO'])+1
    
    #diversity = rec.nucl_diversity
    #if diversity is None:
    #    diversity=-1
    diversity=-1    

    fvals=list()
    for f in filters:
        if f in rec.FILTER:
            fvals+=[1]
        else:
            fvals+=[0]

    DD2MAF=-1
    if 'DD2FDKMAF' in rec.INFO:
        DD2MAF=rec.INFO['DD2FDKMAF']



    print >>sys.stdout, "\t".join(map(str,[
                    rec.CHROM,
                    rec.POS,
                    vartype,
#                    '.'.join(alleles),
                    varlen,
                    bcomplex,
                    STRtype,
                    INDtype,
                    dust,
                    STRP,
                    STRE,
                    alleleCount,
                    round(rec.call_rate,3),
                    round(diversity,3)]+
                    fvals+
                    [coding,
                    consequence,
                    DD2MAF,
                    ]))
                     

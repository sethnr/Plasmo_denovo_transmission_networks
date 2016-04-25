#!/usr/bin/python

import vcf
import sys
import argparse
from string import *
import subprocess
from os.path import basename, splitext
from collections import Counter

parser = argparse.ArgumentParser(description='get STR length differences from disco VCF')

parser.add_argument('-v','--vcf', dest='vcfs', action='append', help='all vcf files for comparison')

parser.add_argument('-o','--out', action="store", dest='outFile', type=str, help='make combined outfile named <outfile> instead of separate files', nargs='?', default=None)
#parser.add_argument('-i','--defaultSize', action="store", dest='defSize', type=int, help='default size for STRs', nargs='?', default=0)


args = parser.parse_args()

vcfs = []
vcfFiles = []
vcfNames = []
vcfSamples = []
sampN = 0

print args.vcfs

vcfSNPs = dict()
vcfINDELs = dict()
allSNPs=set()
allINDELs=set()
for vcfname in args.vcfs:
    print vcfname
    vcfindels =set()
    vcfsnps=set()
    vcffile = open(vcfname,'r')
    reader=vcf.Reader(vcffile)
    i=0
    print reader
    for rec in reader:
        alls = '\t'.join(map(str,[rec.REF]+rec.ALT))
#        print "\t".join(map(str,[rec.CHROM,rec.POS,alls]))
        if rec.is_indel:
            vartype="INDEL"
            vcfindels.add((rec.CHROM,rec.POS,alls))
        elif rec.is_snp:
            vartype = "SNP"
            vcfsnps.add((rec.CHROM,rec.POS,alls))
    vcfSNPs[vcfname]=vcfsnps
    vcfINDELs[vcfname]=vcfindels
    allSNPs = allSNPs.union(vcfsnps)
    allINDELs = allINDELs.union(vcfindels)



venn=dict()    
for snp in allSNPs:
    (chrom,pos,alls) = snp
    print "\t".join([chrom,str(pos),alls,"SNP"]),
    combs = ""
    for vcffile in args.vcfs:
        inset = '0'
        if snp in vcfSNPs[vcffile]:
            inset='1'
        print "\t"+inset,
        combs+=inset
    if ('SNP',combs) not in venn: venn[('SNP',combs)]=0
    venn[('SNP',combs)]+=1
    print ""
    
for indel in allINDELs:
    (chrom,pos,alls) = indel
    print "\t".join([chrom,str(pos),alls,"INDEL"]),
    combs=""
    for vcffile in args.vcfs:
        inset='0'
        if indel in vcfINDELs[vcffile]:
            inset='1'
        print "\t"+inset, 
        combs+=inset
    if ('INDEL',combs) not in venn: venn[('INDEL',combs)]=0
    venn[('INDEL',combs)]+=1
    print ""


for vtype,comb in sorted(venn):
    print >>sys.stderr, vtype,comb,venn[(vtype,comb)]

#!/usr/bin/python

import sys
import argparse
from string import *
from subprocess import call, check_output
import os.path as path
import os
from math import ceil,floor
import pysam
import vcf
import re
import primer3 as p3
from pyfaidx import Fasta
import tabix 
import copy

primvals={
    'PRIMER_TASK' : 'generic',
#    'PRIMER_PICK_INTERNAL_OLIGO':0,
    'PRIMER_PICK_LEFT_PRIMER':1,
    'PRIMER_PICK_RIGHT_PRIMER':1,
#    'PRIMER_PICK_ANYWAY': 1,
    'PRIMER_LOWERCASE_MASKING': 1,
    'PRIMER_OPT_SIZE': 20,
    'PRIMER_MIN_SIZE': 18,
    'PRIMER_MAX_SIZE': 25,
    'PRIMER_OPT_TM': 60.0,
    'PRIMER_MIN_TM': 57.0,
    'PRIMER_MAX_TM': 63.0,
    'PRIMER_MIN_GC': 20.0,
    'PRIMER_MAX_GC': 80.0,
    'PRIMER_MAX_POLY_X': 100,
    'PRIMER_SALT_MONOVALENT': 50.0,
    'PRIMER_DNA_CONC': 50.0,
    'PRIMER_MAX_NS_ACCEPTED': 0,
    'PRIMER_MAX_SELF_ANY': 12,
    'PRIMER_MAX_SELF_END': 8,
    'PRIMER_PAIR_MAX_COMPL_ANY': 12,
    'PRIMER_PAIR_MAX_COMPL_END': 8,
    'PRIMER_EXPLAIN_FLAG': 1,
    'PRIMER_PRODUCT_OPT_SIZE':50,
    'PRIMER_PRODUCT_SIZE_RANGE': [[50,100]]
}



parser = argparse.ArgumentParser(description='identify primer pairs for STR aplification')

parser.add_argument('-v','--vcf', action="store", dest='vcf', type=str, help='vcf file of INDEL calls', nargs='?', default=None)
parser.add_argument('-f','--fasta', action="store", dest='fasta', type=str, help='fasta file (ref)', nargs='?', default=None)
parser.add_argument('--fout','--fasta_out', action="store", dest='fastaOut', type=str, help='fasta file for dusted sequence', nargs='?', default=None)
parser.add_argument('-d','--dust', action="store", dest='dust', type=str, help='dust file (1/0 file of low-complexity regions)', nargs='?', default=None)

parser.add_argument('-o','--out', action="store", dest='outFile', type=str, help='outFile', nargs='?', default=None)


parser.add_argument('-F','--flank', action="store", dest='flank', type=int, help='flank to parse out around refs (bp)', nargs='?', default=1000)
parser.add_argument('-t','--tmp', action="store", dest='tmp', type=str, help='temp folder for calculations <tmp>', nargs='?', default='tmp')


args = parser.parse_args()

TMP="./"+args.tmp+"/"
if not path.exists(TMP):
    os.makedirs(TMP)

name = path.basename(args.fasta).replace('.fasta','')
dustfile = tabix.open(args.dust)
if args.fastaOut is not None:
    fastaOut = open(args.fastaOut,'w')
fasta = Fasta(args.fasta)

def _parse_chrlens_from_dict(fasta):
    seqdictF = fasta.replace(".fasta",".dict")
    seqdict = open(seqdictF,'r')
    locs = []
    lengths = dict()
    for line in seqdict:
        snres = re.search('SN:(\S+)', line)
        if snres is not None:
            seqname = snres.group(1)
            seqlen = re.search('LN:(\d+)', line).group(1)
#            print >>sys.stderr, seqname, seqlen
            lengths[seqname] = int(seqlen)
    return lengths

def _getDustRegion(chrom, st, en):
#    print >>sys.stderr, "reading in dust file "+args.dust+" : "+chrom
    dust = dict()
    for line in dustfile.query(chrom,st,en):
        (chrom, start, end) = line
        start = int(start)
        end = int(end)
        if end > en: end =en
        for p in range(start,end):
            dust[(chrom,p)] = (chrom,start, end)
#    print >>sys.stderr, str(len(dust))+" bases dusted in "+chrom+":"+str(st)+"-"+str(en)
    
    return sorted(dust)

def _merge_regions(locs, joinflank=0):
    i=1
    locs = sorted(locs)
    while i < len(locs):
        #print i
        (c1, s1, e1) = locs[i-1]
        (c2, s2, e2) = locs[i]
        if (c1 == c2) and (e1 + joinflank > s2):
        #    print "removing",locs[i-1],locs[i]
            locs.remove(locs[i])
            locs.remove(locs[i-1])
            locs += [(c1, s1, e2)]
            locs = sorted(locs)
        else: 
            i += 1
    return locs


def _get_dusted_seq(c,st,en):
    dust =_getDustRegion(c,st,en)    
    seq = list(fasta[c][st:en].seq)
    smseq = list(fasta[c][st:en].seq)
#    print >>sys.stderr, len(dust),len(seq)
    for c,p in dust:
#        print c,p
        i = p-st-1
        if i > 0 and i <= len(seq):
            seq[i]='n'
            smseq[i] = smseq[i].lower()
    print >>sys.stderr, str(round(float(len(dust))/len(seq),2))+"% flanking seq dusted"
    smseq = str(''.join(smseq))
    print >>sys.stderr, smseq
    seq = str(''.join(seq))
    return smseq


def _find_nearest_primers(var, flank=1000, increment=100):
    vst = var["vstart"]
    ven=var["vend"]
    vrlen=var["vend"]-var["vstart"]
    c =var["chrom"]
    regst = vst-flank
    regen = ven+flank-1
    #get dusted sequence
    name=c+":"+str(vst)
#    print >>sys.stderr, c, regst, regen
    if regst <= 0: 
        flank=flank+regst
        regst=1
    if regen > chrlens[c]: 
        regen=chrlens[c]-1
#    print >>sys.stderr, c, regst, regen

    seq = _get_dusted_seq(c,regst,regen)
#    print >>sys.stderr, len(seq)

    if args.fastaOut is not None:
        print >>fastaOut, "\t".join(map(str,[name,regst,regen]))
        print >>fastaOut, seq
        
    seqvals={'SEQUENCE_ID': name,
             'SEQUENCE_TEMPLATE': seq,
#             'SEQUENCE_INCLUDED_REGION': [1,(2*flank)],
             'SEQUENCE_INCLUDED_REGION': [0,len(seq)],
             'SEQUENCE_TARGET': [flank,vrlen]
}
    #print >>sys.stderr, seqvals

    p3res=dict()

    i="0"

    myprimvals=copy.copy(primvals)
    while not 'PRIMER_PAIR_'+i+'_PRODUCT_SIZE' in p3res:
        p3res = p3.bindings.designPrimers(seqvals,myprimvals)
#        sizeRange = primvals['PRIMER_PRODUCT_SIZE_RANGE'][0]
#        print >>sys.stderr, sizeRange
        minS,maxS = myprimvals['PRIMER_PRODUCT_SIZE_RANGE'][0]
        myprimvals['PRIMER_PRODUCT_SIZE_RANGE'] = [[maxS,maxS+increment]]
        increment += increment
        print >>sys.stderr, maxS
        if minS > flank: break
        if maxS > len(seq): break
                                                 
    #for i in range(0,5):
    #    i = str(i)
    #    print >>sys.stderr, i
#    for key in sorted(p3res):
#        print key, p3res[key]
    if not 'PRIMER_PAIR_'+i+'_PRODUCT_SIZE' in p3res:
        print >>sys.stderr, "no primer found within "+str(flank)+"bp"
#        print >>sys.stderr, p3res["PRIMER_PAIR_EXPLAIN"]
        return -1,-1,-1,"","",-1,-1
    else:
        pLpos, pLlen = map(int,p3res['PRIMER_LEFT_'+i])
        pRpos, pRlen = map(int,p3res['PRIMER_RIGHT_'+i])
        pSize = p3res['PRIMER_PAIR_'+i+'_PRODUCT_SIZE']
        pLpos = pLpos + regst
        pRpos = pRpos + regst
        pL = p3res['PRIMER_LEFT_'+i+"_SEQUENCE"]
        pR = p3res['PRIMER_RIGHT_'+i+"_SEQUENCE"]
#        print >>sys.stderr, p3res
        pLtm = p3res['PRIMER_LEFT_'+i+'_TM']
        pRtm = p3res['PRIMER_RIGHT_'+i+'_TM']
        return pLpos, pRpos, pSize, pL, pR, pLtm, pRtm
#        print >>sys.stderr, pLpos, pLlen, pSize
#        print >>sys.stderr, pRpos, pRlen, pSize
    
    
    
    


chrlens =  _parse_chrlens_from_dict(args.fasta)

vcfFile = open(args.vcf,'r')
reader=vcf.Reader(vcfFile)

if args.outFile is not None:
    outfile = open(args.outFile,'w')
else:
    outfile = sys.stdout
rvars=dict()
rfilts=dict()
for rec in reader:
    #print rec.CHROM,str(rec.POS),str(rec.is_indel)
    if not rec.is_indel:
        continue
    var = dict()
    rlen = len(rec.REF)
    vstart = rec.POS
    vend = rec.POS + rlen
    var["chrom"] = rec.CHROM
    var["vstart"] = vstart
    var["vend"] = vend
    var["rlen"] = rlen
    var["vlen"] = (rlen - len(a) for a in rec.ALT)
    
    pchrom = rec.CHROM
    pstart = vstart-args.flank
    pend = vend+args.flank
    #regions += [(pchrom,pstart,pend)]

    pStart,pEnd,prodSize,p1,p2 = -1,-1,-1,"",""

    srchFlank = 0
    maxFlank=args.flank
    while pStart == -1 and srchFlank < maxFlank:
        srchFlank += args.flank
        print >>sys.stderr, srchFlank
        pStart,pEnd,prodSize,p1,p2,pt1,pt2 = _find_nearest_primers(var,srchFlank)
#    if prodSize > -1:
#        pSizes = [prodSize] + [prodSize+l for l in var["vlen"]]
#        pSizes = map(str,pSizes)
#    else: 
#        pSizes = ["-1"]
    
    reflen = len(rec.REF)
    vSizes=[]
    for call in rec.samples:
        bases = call.gt_bases
#        print >>sys.stderr, bases
        if bases is not None:
            GT = bases.split("/")[0]
            vsize = len(GT)-reflen
        else:
            vsize= -1
        vSizes += [vsize]

    pSizes = [v+prodSize for v in vSizes]
    
    var["pStart"]=pStart
    var["pEnd"]=pEnd
    var["prodSize"]=prodSize
    var["vSizes"]=vSizes
    var["prodSizes"]=pSizes
    var["p5"]=p1
    var["p3"]=p2
    var["pt5"]=pt1
    var["pt3"]=pt2
    if (var["chrom"],pStart,pEnd) not in rvars:
        rvars[(var["chrom"],pStart,pEnd)]=[]
    rvars[(var["chrom"],pStart,pEnd)] += [var]

    if (var["chrom"],pStart,pEnd) not in rfilts:
        rfilts[(var["chrom"],pStart,pEnd)]=[]
    for f in rec.FILTER:
        rfilts[(var["chrom"],pStart,pEnd)] += [f]

#    print >>outfile, "\t".join(map(str,[var["chrom"],
    print >>sys.stderr, "\t".join(map(str,[var["chrom"],
                    var["vstart"],
                    #rec.REF,
                    pStart,
                    pEnd,
                    "/".join(map(str,pSizes)),
                    #prodSize,
                    #prodSize+max(var["vlen"]),
                    p1,p2])
                    )



#regions = _merge_regions(regions,mergeFlank=args.flank)
#
for c,s,e in sorted(rvars):
    varlist = rvars[(c,s,e)]
    outvar = varlist[0]
    a = outvar["vSizes"]
    outvar["vstart"] = [outvar["vstart"]]
    for var in varlist[1:]:
        b = var["vSizes"]
        outvar["vSizes"] = [a[i]+b[i] for i in range(0,len(a))]
        outvar["vstart"] += [var["vstart"]]
    outvar["prodSizes"] = [outvar["prodSize"]+v for v in [0]+outvar["vSizes"]]
    pcDiff = (round((max(outvar["prodSizes"])/float(min(outvar["prodSizes"]))),3)-1)

    print >>outfile, "\t".join(map(str,[outvar["chrom"],
                    "/".join(map(str,outvar["vstart"])),
                    outvar["pStart"],
                    outvar["pEnd"],
                    ",".join(list(set(rfilts[(c,s,e)]))),
                    #"/".join(map(str,outvar["vSizes"])),
                    "/".join(map(str,outvar["prodSizes"])),
                    pcDiff,
                    outvar["p5"],outvar["p3"],
                    outvar["pt5"],outvar["pt3"]])
                    )

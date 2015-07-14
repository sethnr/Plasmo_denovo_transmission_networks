#!/usr/bin/python

import vcf
import sys
import argparse
from string import *

parser = argparse.ArgumentParser(description='get STR length differences from disco VCF')

parser.add_argument('--v1','--vcf1', action="store", dest='vcfFile1', type=str, help='vcfFile1', nargs='?', default=None)
parser.add_argument('--v2','--vcf2', action="store", dest='vcfFile2', type=str, help='vcfFile2', nargs='?', default=None)
parser.add_argument('-n','--nucmer', action="store", dest='nucmer', type=str, help='nucmer', nargs='?', default=None)

parser.add_argument('-o','--out', action="store", dest='outFile', type=str, help='outFile', nargs='?', default=None)
#parser.add_argument('-i','--defaultSize', action="store", dest='defSize', type=int, help='default size for STRs', nargs='?', default=0)


args = parser.parse_args()
print >>sys.stderr, args.vcfFile1, args.vcfFile2


outfile = sys.stdout
if args.outFile is not None:
    outfile = open(args.outFile,'w')

    
vcfFile1 = open(args.vcfFile1,'r')
reader1=vcf.Reader(vcfFile1)

vcfFile2 = open(args.vcfFile2,'r')
reader2=vcf.Reader(vcfFile2)


vcfoutF1 = replace(args.vcfFile1,'.vcf','.cfout.vcf')
vcfoutF1 = replace(vcfoutF1,'.vcf.gz','.cfout.vcf')
print >>sys.stderr, "VCFOUT1",vcfoutF1
vcfoutF1 = open(vcfoutF1,'w')
vcfout1=vcf.Writer(vcfoutF1,reader1)

vcfoutF2 = replace(args.vcfFile2,'.vcf','.cfout.vcf')
vcfoutF2 = replace(vcfoutF2,'.vcf.gz','.cfout.vcf')
print >>sys.stderr, "VCFOUT2",vcfoutF2
vcfoutF2 = open(vcfoutF2,'w')
vcfout2=vcf.Writer(vcfoutF2,reader2)


v1end = False
v2end = False

c1=None  #chr vcf1
c2=None  #chr vcf2
p1=-1    #pos vcf1
p2=-1    #pos vcf2
ar1=None #ref allele vcf1
aa1=None #alt allele vcf1
ar2=None #ref allele vcf2
aa2=None #alt allele vcf2
calls1=[] #calls vcf1
calls2=[] #calls vcf2

def _readVar(reader):
    try:
        rec = reader.next()
    except StopIteration:
        return(None, None, -1, None, None, [], None, True) 
    calls = []
    type = "SNP"
    if len(rec.REF) > 1 or len(rec.ALT) > 1:
        type="INDEL"
    for call in rec.samples:
        calls += [call.gt_type]
    return (rec, rec.CHROM, rec.POS, rec.REF, rec.ALT, calls, type, False)
    
sameRef = True

pv1=0
pv2=0
pi1=0
ps1=0
pi2=0
ps2=0
match=0
matchI=0
matchS=0
mismatch=0


def _handlePriv1():
    global pv1, ps1, pi1, rec1
    pv1 +=1
    if type1 == "SNP":
        ps1 +=1
    elif type1 == "INDEL":
        pi1 +=1
    rec1.INFO['TYPE']=type1
    rec1.INFO['CON']='PRIVATE'
    vcfout1.write_record(rec1)

def _handlePriv2():
    global pv2, ps2, pi2, rec2
    pv2 +=1
    if type2 == "SNP":
        ps2 +=1
    elif type2 == "INDEL":
        pi2 +=1
    rec2.INFO['TYPE']=type2
    rec2.INFO['CON']='PRIVATE'
    vcfout2.write_record(rec2)

def _handleMatch():
    global matchI,match,matchS, type1, type2, rec1, rec2
    
    if type2 == "SNP":
        matchS +=1
    elif type2 == "INDEL":
        matchI +=1
    rec2.INFO['TYPE']=type2
    rec2.INFO['CON']='MATCH'
    rec1.INFO['TYPE']=type1
    rec1.INFO['CON']='MATCH'
    vcfout1.write_record(rec1)
    vcfout2.write_record(rec2)

def _handleMismatch():
    global mismatch, type1, type2, rec1, rec2
    rec2.INFO['TYPE']=type2
    rec2.INFO['CON']='MISMATCH'
    rec1.INFO['TYPE']=type1
    rec1.INFO['CON']='MISMATCH'
    vcfout1.write_record(rec1)
    vcfout2.write_record(rec2)



if sameRef:
    while not (v1end and v2end):
        #get initial values
        print >> sys.stderr, c1, p1, c2, p2,
        if c1 is None and c2 is None:
            #neither file has started
            rec1, c1, p1, ar1, aa1, calls1, type1, v1end = _readVar(reader1)
            rec2, c2, p2, ar2, aa2, calls2, type2, v2end = _readVar(reader2)
        elif v2end:
            #file 2 has finished
            print >>sys.stderr, '<'
##             pv1 +=1
##             if type1 == "SNP":
##                 ps1 +=1
##             elif type1 == "INDEL":
##                 pi1 +=1
            _handlePriv1()
            # get next
            rec1, c1, p1, ar1, aa1, calls1, type1, v1end = _readVar(reader1)
        elif v1end:
            #file 1 has finished
            print >>sys.stderr, '>'
##             pv2 +=1
##             if type2 == "SNP":
##                 ps2 +=1
##             elif type2 == "INDEL":
##                 pi2 +=1
            _handlePriv2()
            rec2, c2, p2, ar2, aa2, calls2, type2, v2end = _readVar(reader2)
            
        elif (c1 == c2 and p1 > p2) or c1 > c2:
            # print cf
            print >>sys.stderr,'>'
##             pv2 +=1
##             if type2 == "SNP":
##                 ps2 +=1
##             elif type2 == "INDEL":
##                 pi2 +=1
            _handlePriv2()
            # get next
            rec2, c2, p2, ar2, aa2, calls2, type2, v2end = _readVar(reader2)
        elif (c1 == c2 and p1 < p2) or c1 < c2:
            # print cf
            print >>sys.stderr,'<'
##             pv1 +=1
##             if type1 == "SNP":
##                 ps1 +=1
##             elif type1 == "INDEL":
##                 pi1 +=1
            _handlePriv1()
            # get next
            rec1, c1, p1, ar1, aa1, calls1, type1, v1end = _readVar(reader1)
        elif p1 == p2 and c1 == c2:
            # print cf
            if ar1 == ar2 and aa2 == aa1:
                #match = True
                print >>sys.stderr,"="
##                 match +=1
##                 if type2 == "SNP":
##                     matchS +=1
##                 elif type2 == "INDEL":
##                     matchI +=1
                _handleMatch()
            else:
##                mismatch += 1
                print >>sys.stderr,"."
                _handleMismatch()
            # get next
            rec1, c1, p1, ar1, aa1, calls1, type1, v1end = _readVar(reader1)
            rec2, c2, p2, ar2, aa2, calls2, type2, v2end = _readVar(reader2)
        else:
            # if none have worked, we got issues (end of file)
            print >>sys.stderr,"!"

elif recipRef:
    pass
##         elif p1 == pTrans(p2) and c1 == cTrans(c2):
##             if ar1 == aa2 and ar2 == aa1:
##                 match = True
##                 print >>sys.stderr,"="
##             else:
##                 print >>sys.stderr,"."
            # get next


print "pv1", "pv2", "match", "mismatch"
print "vars:",pv1, pv2, match, mismatch
print "snps:",ps1, ps2, matchS
print "indels:",pi1, pi2, matchI

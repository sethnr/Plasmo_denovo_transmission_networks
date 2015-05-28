#!/usr/bin/python

import sys
from Bio import SeqIO
from Bio.Seq import Seq
import os.path as path
import argparse
import vcf
from subprocess import call, check_output


########
# get some vars
########
parser = argparse.ArgumentParser(description='build index of genome matches in 3k blocks')
parser.add_argument('-f','--fasta', action="store", dest='fasta', type=str, help='target genome to index (fasta)', nargs='?')
parser.add_argument('-m','--mummer', action="store", dest='mummer', type=str, help='bed file for SNP locations', nargs='?', default=None)
parser.add_argument('--v1','--vcf1', action="store", dest='vcfFile1', type=str, help='vcfFile', nargs='?', default=None)
parser.add_argument('--v2','--vcf2', action="store", dest='vcfFile2', type=str, help='vcfFile', nargs='?', default=None)
parser.add_argument('-v','--vcf', action="store", dest='vcf', type=str, help='vcf template (for metadata)', nargs='?', default=None)
parser.add_argument('--sanity', action="store_true", dest='sanity', help='print actual seq from file to check for insert seq', default=False)

args = parser.parse_args()


#get coordinates (start/end/seq) for each ref
global seq1, seq2, st1, st2, en1, en2, lc1, lc2, inseq, delseq

seq1='.'
seq2='.'
st1=0
st2=0
en1=0
en2=0

lc1=None
lc2=None
inseq=''
delseq=''


#defaults
qual=40
filter='PASS'
info='.'
format="GT"

def _write_indel():
    global seq1, seq2, st1, st2, en1, en2, lc1, lc2, inseq, delseq
    print >>sys.stderr, i," ", lc1, st1, en1, inseq, lc2, st2, en2, delseq
    #check if indel is in memory
    if inseq !='' or delseq != '':
    
        pb1 = _get_prevbase(fasta1,lc1,st1)
        pb2 = _get_prevbase(fasta2,lc2,st2)

        if inseq != '':
            delseq=pb1+delseq
            inseq=pb2+inseq
            print >>sys.stderr, "  IN ",lc1, st1, en1, inseq, lc2, st1, en2, '.'

            # pyvcf not great for writing new vcfs
            # try a different approach...
            #line = vcf.model._Record(lc1, st1, '.', delseq, inseq, qual, filter, info, format,[1,2], [0,1])# [vcfname1,vcfname2])
            #vcf1.write_record(line)
            vcf1line = [lc1, st1, '.', delseq, inseq, qual, filter, "TYPE=ins", format, 0, 1]
            print >>vcf1, "\t".join(map(str,vcf1line))
            vcf2line = [lc2, st2, '.', inseq, delseq, qual, filter, "TYPE=del", format, 0, 1]
            print >>vcf2, "\t".join(map(str,vcf2line))
            
        elif delseq != '':
            delseq=pb1+delseq
            inseq=pb2+inseq
            print >>sys.stderr, "  DEL", lc1, st1, en1, '.',   lc2, st1, en2, delseq
            vcf1line = [lc1, st1, '.', delseq, inseq, qual, filter, "TYPE=del", format, 0, 1]
            print >>vcf1, "\t".join(map(str,vcf1line))
            vcf2line = [lc2, st2, '.', inseq, delseq, qual, filter, "TYPE=ins", format, 0, 1]
            print >>vcf2, "\t".join(map(str,vcf2line))

#        if args.sanity: print >>vcf1, "\t".join(map(str,[lc1,st1,en1,fasta1,getseq(fasta1,lc1,st1,en1)]))
#        if args.sanity: print >>vcf2, "\t".join(map(str,[lc2,st2,en2,fasta2,getseq(fasta2,lc2,st2,en2)]))
        if args.sanity: print >>vcf1, "\t".join(map(str,[lc1,st1,en1,fasta1]))
        if args.sanity: print >>vcf2, "\t".join(map(str,[lc2,st2,en2,fasta2]))
 
    #pass

#def _write_del():
    
    #pass

def _get_prevbase(fasta, chrom, pos):
#    if pos2 is None: pos2=pos
#    command = ['samtools', 'faidx', fasta, chrom+":"+str(pos-1)+"-"+str(pos2)]
#    prev = check_output(command)
    #print >>sys.stderr, prev
#    prevbase = prev.splitlines()[1]
    #print prevbase+"!"
    prevbase = getseq(fasta, chrom, pos-1)
    return prevbase
#    faidx = check_output(command)
    
def getseq(fasta, chrom, pos, pos2=None):
    if pos2 is None: pos2=pos
    command = ['samtools', 'faidx', fasta, chrom+":"+str(pos)+"-"+str(pos2)]
    prev = check_output(command)
    #print >>sys.stderr, prev
    prevbase = prev.splitlines()[1]
    #print prevbase+"!"
    return prevbase
#    faidx = check_output(command)
    
def _reset_indels():
    #print >>sys.stderr, "resetting"
    global seq1, seq2, st1, st2, en1, en2, lc1, lc2, inseq, delseq
    seq1='.'
    seq2='.'
    st1=-1
    st2=-1
    en1=-1
    en2=-1

    lc1=None
    lc2=None
    inseq=''
    delseq=''

def _set_starts(c1, c2, s1, s2):
    global st1, st2, lc1, lc2
    st1=s1
    st2=s2
    lc1=c1
    lc2=c2
    

mummer = open(args.mummer)
vcf1 = None
vcf2 = None

#vcfHeader = vcf.Reader(filename=args.vcf)

i=0
for mumm in mummer:
    i+=1
    # [P1]	[SUB]	[SUB]	[P2]	[BUFF]	[DIST]	[FRM]	[TAGS]
    # buff = distance to nearest other SNP
    # dist = dist to end of ref
    # FRM = alignment direction 1=forward, -1=ref
#    print len(mumm.split())

    mumm = mumm.split()

    if(len(mumm)==2):
        fasta1, fasta2 = mumm
        if args.vcfFile1 is None:
            vcfname1 = path.splitext(path.basename(fasta1))[0]
        else:
            vcfname1= args.vcfName1

        if args.vcfFile2 is None:
            vcfname2 = path.splitext(path.basename(fasta2))[0]
        else:
            vcfname2=args.vcfName2
        print >>sys.stderr, "opening vcfs:",vcfname1,"/", vcfname2
        #do not use pyvcf, not great at making new VCFs
        #vcf1 = vcf.Writer(open(vcfname1+".vcf", 'w'),vcfHeader)
        #vcf2 = vcf.Writer(open(vcfname2+".vcf", 'w'),vcfHeader)
        vcf1 = open(vcfname1+".vcf", 'w')
        vcf2 = open(vcfname2+".vcf", 'w')

        
    elif(len(mumm)<10):
        pass
    else:  
        (p1, s1, s2, p2, buff, dist, fr, this, c1, c2) = mumm
        p1 = int(p1); p2 = int(p2)

        #if something is in buffer
        if lc1 is None:
            _set_starts(c1,c2,p1,p2)
            
        elif s1 != '.' and s2 != '.':
            #is a SNP
            #if inseq !='' or delseq != '':
            #if holding indel, write indel

            print >>sys.stderr, "in SNP, writing last indel",
            _write_indel()
            _reset_indels()
            _set_starts(c1,c2,p1,p2)
             
            #pass
        else:
            # is an indel
            # if indel, but coords not consecutuve
            if inseq != '' and ((en2+1) != p2 or lc1 != c1):
                print >>sys.stderr,"insert break (",en1,"+1 != ",p1,")",
                #record insertion
                _write_indel()
                #reset values
                _reset_indels()
                #put in starter vals
                _set_starts(c1,c2,p1,p2)

            # if deletion, but coords not consecutuve
            elif delseq != '' and ((en1+1) != p1 or lc1 != c1):
                print >>sys.stderr, "del break (",en1,"+1 != ",p1,")",
                #record deletion
                _write_indel()
                #reset values
                _reset_indels()
                #put in starter vals
                _set_starts(c1,c2,p1,p2)
            else:
                #continuing indel, just increment previous lines
                pass
             
        #        print "continuing"
        #print s1, s2, inseq, delseq
        #assess current:
        lc1=c1
        lc2=c2
        en1=int(p1)
        en2=int(p2)
        if s2=='.':
            #deletion
            delseq += s1
        #    e2 = p2
        if s1=='.':
            inseq += s2
        #    e1 = p1

i=-1
_write_indel()

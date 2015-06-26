#!/usr/bin/python

import sys
from Bio import SeqIO
from Bio.Seq import Seq
import os.path as path
import argparse
import vcf
from subprocess import call, check_output
import re

########
# get some vars
########
parser = argparse.ArgumentParser(description='build index of genome matches in 3k blocks')
parser.add_argument('-f','--fasta', action="store", dest='fasta', type=str, help='target genome to index (fasta)', nargs='?')
parser.add_argument('-m','--mummer', action="store", dest='mummer', type=str, help='bed file for SNP locations', nargs='?', default=None)
parser.add_argument('--v1','--vcf1', action="store", dest='vcfFile1', type=str, help='vcfFile', nargs='?', default=None)
parser.add_argument('--v2','--vcf2', action="store", dest='vcfFile2', type=str, help='vcfFile', nargs='?', default=None)
parser.add_argument('-v','--vcf', action="store", dest='vcf', type=str, help='vcf template (for metadata)', nargs='?', default='')
parser.add_argument('--sanity', action="store_true", dest='sanity', help='print actual seq from file to check for insert seq', default=False)

args = parser.parse_args()


#get coordinates (start/end/seq) for each ref
global seq1, seq2, st1, st2, en1, en2, lc1, lc2 #, inseq, delseq

seq1=''
seq2=''
st1=0
st2=0
en1=0
en2=0
nearestSNP=-1

lc1=None
lc2=None
#inseq=''
#delseq=''


vartype=''

#defaults
qual=40
filter='PASS'
info='.'
format="GT"

def _write_var():
    global seq1, seq2, st1, st2, en1, en2, lc1, lc2, vartype, lastFR#, inseq, delseq
#    print >>sys.stderr, i," ", lc1, st1, en1, inseq, lc2, st2, en2, delseq
#    print >>sys.stderr, i," ", lc1, st1, en1, seq1, lc2, st2, en2, seq2, vartype, lastFR
    #check if indel is in memory
#    if inseq !='' or delseq != '':

    var2type='SNP'
    pb1=''
    pb2=''
    
    if vartype == 'SNP':
        seq1_1=seq1
        seq2_1=seq2
#        print >>sys.stderr, "VARTYPE = SNP, lastFR=",lastFR
        if lastFR==1:
#            print >>sys.stderr, "FORWARD"
            seq1_2=seq1
            seq2_2=seq2
        elif lastFR==-1:
#            print >>sys.stderr, "REVERSE"
            seq1_2=revcomp(seq1)
            seq2_2=revcomp(seq2)
        
    else:
        pb1 = _get_base(fasta1,lc1,st1-1)
        pb2 = _get_base(fasta2,lc2,st2-1)
        nb1 = _get_base(fasta1,lc1,en1+1)
        nb2 = _get_base(fasta2,lc2,en2+2)
        st2 = st2-1 # count both from 
        st1 = st1-1 # previous base
        
        
        rev=False

        #SEQ 1 and 2, as in first VCF
        seq1_1 = pb1+seq1
        seq2_1 = pb1+seq2

        if vartype=='ins' and lastFR==1:
            var2type='del'
            en1 = en1-1
            #SEQ 1 and 2, as in second VCF
            seq1_2 = pb2+seq1
            seq2_2 = pb2+seq2

        if vartype=='ins' and lastFR==-1:
            var2type='del'
            en1=en1-1
            rev=True
            seq1_2 = pb2+revcomp(seq1)
            seq2_2 = pb2+revcomp(seq2)

        elif vartype=='del' and lastFR==1:
            var2type='ins'
            en2=en2-1
            seq1_2 = pb2+seq1
            seq2_2 = pb2+seq2

        elif vartype=='del' and lastFR==-1:
            var2type='ins'
            en2=en2-1
            rev=True
            seq1_2 = pb2+revcomp(seq1)
            seq2_2 = pb2+revcomp(seq2)

        if rev is True:
            vartype += 'R'
            var2type += 'R'
        else:
            vartype += 'F'
            var2type += 'F'


##         if inseq != '':
##             delseq=pb1+delseq
##             inseq=pb2+inseq
##             print >>sys.stderr, "  IN ",lc1, st1, en1, inseq, lc2, st1, en2, '.'

##             # pyvcf not great for writing new vcfs
##             # try a different approach...
##             #line = vcf.model._Record(lc1, st1, '.', delseq, inseq, qual, filter, info, format,[1,2], [0,1])# [vcfname1,vcfname2])
##             #vcf1.write_record(line)

    if args.sanity:  #PUT IN INDEX LINE TO ALLOW REORDERING OF FILE
        print >>vcf1,i,"\t",
        print >>vcf2,i,"\t",
        
    
    vcf1line = [lc1, st1, '.', seq1_1, seq2_1, qual, filter, "TYPE="+vartype, format, 0, 1]
    print >>vcf1, "\t".join(map(str,vcf1line))
    vcf2line = [lc2, st2, '.', seq2_2, seq1_2, qual, filter, "TYPE="+var2type, format, 0, 1]
    print >>vcf2, "\t".join(map(str,vcf2line))
            
##         elif delseq != '':
##             delseq=pb1+delseq
##             inseq=pb2+inseq
##             print >>sys.stderr, "  DEL", lc1, st1, en1, '.',   lc2, st1, en2, delseq
##             vcf1line = [lc1, st1, '.', delseq, inseq, qual, filter, "TYPE=del", format, 0, 1]
##             print >>vcf1, "\t".join(map(str,vcf1line))
##             vcf2line = [lc2, st2, '.', inseq, delseq, qual, filter, "TYPE=ins", format, 0, 1]
##             print >>vcf2, "\t".join(map(str,vcf2line))
    
    fsize=20
    varflank1_1 = getseq(fasta1,lc1,st1-fsize,st1-1).lower()+seq1_1+getseq(fasta1,lc1,en1+1,en1+fsize).lower()
    varflank2_1 = getseq(fasta2,lc2,st2-fsize,st2-1).lower()+seq2_1+getseq(fasta2,lc2,en2+1,en2+fsize).lower()
    varflank1_2 = getseq(fasta1,lc1,st1-fsize,st1-1).lower()+seq1_2+getseq(fasta1,lc1,en1+1,en1+fsize).lower()
    varflank2_2 = getseq(fasta2,lc2,st2-fsize,st2-1).lower()+seq2_2+getseq(fasta2,lc2,en2+1,en2+fsize).lower()

    if args.sanity:
        print >>vcf1, "\t".join(map(str,[i, lc1,st1,en1,st2,en2,vartype,buff,"vcf",
                                         varflank1_1,
                                         varflank2_1,
                                         ]))
        print >>vcf1, "\t".join(map(str,[i, lc1,st1,en1,st2,en2,vartype,buff,"seq",
                                         getseq(fasta1,lc1,st1-fsize,en1+fsize).lower(),
                                         #getseq(fasta2,lc2,st2-fsize,en2+fsize).lower(),
                                         ]))
 #       print >>vcf1, "\t".join(map(str,[lc1,st1,en1,st2,en2,vartype,buff,"vcfrev",
 #                                        revcomp(varflank2_1),
 #                                        revcomp(varflank1_1),
 #                                        ]))
        print >>vcf1,""
        
        print >>vcf2, "\t".join(map(str,[i, lc2,st2,en2,st1,en1,var2type,buff,"vfl",
                                         varflank2_2,
                                         varflank1_2,
                                         ]))
        print >>vcf2, "\t".join(map(str,[i, lc2,st2,en2,st1,en1,var2type,buff,"seq",
                                         getseq(fasta2,lc2,st2-fsize,en2+fsize).lower(),
                                         #getseq(fasta1,lc1,st1-fsize,en1+fsize).lower()
                                         ]))
#        print >>vcf2, "\t".join(map(str,[lc2,st2,en2,st1,en1,var2type,buff,"vflrev",
#                                         revcomp(varflank1_1),
#                                         revcomp(varflank2_1),
#                                         ]))
        print >>vcf2,""


    #if args.sanity: print >>vcf1, "\t".join(map(str,[lc1,st1,en1,vartype]))
    #if args.sanity: print >>vcf2, "\t".join(map(str,[lc2,st2,en2,var2type]))
 
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

def _get_base(fasta, chrom, pos):
    prevbase = getseq(fasta, chrom, pos)
    return prevbase

def revcomp(seq):
    d={'a':'t','t':'a','c':'g','g':'c',
       'A':'T','T':'A','C':'G','G':'C'}
    r = [d[b] for b in reversed(seq)]
    return ''.join(r)

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
    global seq1, seq2, st1, st2, en1, en2, lc1, lc2 #, inseq, delseq
    chr1='.'
    chr2='.'
    st1=-1
    st2=-1
    en1=-1
    en2=-1

    lc1=None
    lc2=None
    #inseq=''
    #delseq=''
    seq1=''
    seq2=''
    vartype=''
    nearestSNP=-1

def _set_starts(c1, c2, p1, p2, s1, s2):
    global st1, st2, lc1, lc2, seq1, seq2, vartype
    st1=p1
    st2=p2
    lc1=c1
    lc2=c2
    if   s1!='.' and s2=='.':
        #deletion
        vartype='del'
        seq1 = s1
    elif s1=='.' and s2!='.':
        vartype='ins'
        seq2 = s2
    elif s2!='.' and s2!='.':
        vartype='SNP'
        seq1 = s1
        seq2 = s2

def _setup_fasta(filename, seqdict=None, prefix=''):
    if seqdict==None:
            seqdict = dict()

    samplename = path.splitext(path.basename(filename))[0]
    vcf = open(prefix+samplename+".vcf", 'w')
##     fasta_seqs = check_output(['grep','>',filename])
##     fasta_seqs = fasta_seqs.splitlines()
##     for seq in fasta_seqs:
##         seq = re.search(r">([^\s]*)\s*", seq).group(1)
##         if seq not in seqdict:
##             seqdict[seq]=samplename
##         else:
##             raise KeyError("key "+seq+"("+filename+") already found seq dict")
    return (vcf, samplename, seqdict)
        
mummer = open(args.mummer)

vcf1 = None
vcf2 = None
fasta1 = None
fasta2 = None
seqdict = None
vcfdict = dict()

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

                                               #fasta_name, seqs_in_fasta, prefix for output files
        (vcf1, sample1, seqdict) = _setup_fasta(fasta1, None, args.vcf)
        (vcf2, sample2, seqdict) = _setup_fasta(fasta2, seqdict, args.vcf)
        vcfdict[sample1]=vcf1
        vcfdict[sample2]=vcf2
        
#        print seqdict
##         if args.vcfFile1 is None:
##             vcfname1 = path.splitext(path.basename(fasta1))[0]
##         else:
##             vcfname1= args.vcfName1

##         if args.vcfFile2 is None:
##             vcfname2 = path.splitext(path.basename(fasta2))[0]
##         else:
##             vcfname2=args.vcfName2
##         print >>sys.stderr, "opening vcfs:",vcfname1,"/", vcfname2
##         #do not use pyvcf, not great at making new VCFs
##         #vcf1 = vcf.Writer(open(vcfname1+".vcf", 'w'),vcfHeader)
##         #vcf2 = vcf.Writer(open(vcfname2+".vcf", 'w'),vcfHeader)
##         vcf1 = open(vcfname1+".vcf", 'w')
##         vcf2 = open(vcfname2+".vcf", 'w')

        
    elif(len(mumm)<10):
        pass
    else:  
        (p1, s1, s2, p2, buff, dist, this, fr, c1, c2) = mumm
        p1 = int(p1); p2 = int(p2)

        #if something is in buffer
        if lc1 is None:
            _set_starts(c1,c2,p1,p2,s1,s2)
            
#        elif s1 != '.' and s2 != '.':

        elif lc1 == c1 and lc2 == c2:
            if   ((en1+1)==p1 and (en2==p2)):
                #continue deletion
                seq1 += s1
            elif ((en1==p1)   and (en2+1)==p2):
                #continue insertion
                seq2 += s2
            elif ((en1+1)==p1 and (en2+1)==p2):
                #continue multi line SNP
                seq1 += s1
                seq2 += s2
            else:
                #not a continuation
                #print last var, start new var
                _write_var()
                _reset_indels()
                _set_starts(c1,c2,p1,p2,s1,s2)
           
        else:
            #chromosome break
            #print last var, start new var
            _write_var()
            _reset_indels()
            _set_starts(c1,c2,p1,p2,s1,s2)
           
            

##         elif ((en2+1) != p2 or lc1 != c1):
##             #chromosome break!
##             #is a SNP
##             #if inseq !='' or delseq != '':
##             #if holding indel, write indel

##             print >>sys.stderr, "in SNP, writing last indel",
##             _write_indel()
##             _reset_indels()
##             _set_starts(c1,c2,p1,p2)
             
##             #pass
##         else:
##             # is an indel
##             # if indel, but coords not consecutuve
##             if inseq != '' and ((en2+1) != p2 or lc1 != c1):
##                 print >>sys.stderr,"insert break (",en1,"+1 != ",p1,")",
##                 #record insertion
##                 _write_indel()
##                 #reset values
##                 _reset_indels()
##                 #put in starter vals
##                 _set_starts(c1,c2,p1,p2)

##             # if deletion, but coords not consecutuve
##             elif delseq != '' and ((en1+1) != p1 or lc1 != c1):
##                 print >>sys.stderr, "del break (",en1,"+1 != ",p1,")",
##                 #record deletion
##                 _write_indel()
##                 #reset values
##                 _reset_indels()
##                 #put in starter vals
##                 _set_starts(c1,c2,p1,p2)
##             else:
##                 #continuing indel, just increment previous lines
##                 pass
             
        #        print "continuing"
        #print s1, s2, inseq, delseq
        #assess current:
        lc1=c1
        lc2=c2
        en1=int(p1)
        en2=int(p2)
#        if vartype='':
##         if   s1!='.' and s2=='.':
##             #deletion
##             vartype='del'
##             seq1 += s1
##         elif s1=='.' and s2!='.':
##             vartype='ins'
##             seq2 += s2
##         elif s2!='.' and s2!='.':
##             vartype='SNP'
##             seq1 += s1
##             seq2 += s2
        nearestSNP=buff
        lastFR=int(fr)
i=-1
_write_var()

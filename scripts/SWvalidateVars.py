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


parser = argparse.ArgumentParser(description='SW realign variant-transformed sequence to original reference')

parser.add_argument('-v','--vcf', action="store", dest='vcfFull', type=str, help='vcf file for consensu building', nargs='?', default=None)
parser.add_argument('-r1','-r','--fasta_for_cons', action="store", dest='fasta_cons', type=str, help='fasta file (ref)', nargs='?', default=None)
parser.add_argument('-r2','-f','--fasta_for_align', action="store", dest='fasta_align', help='fasta file (query)', nargs='?')

parser.add_argument('-o','--out', action="store", dest='outFile', type=str, help='outFile', nargs='?', default=None)

#vars to assess
parser.add_argument('-v2','--vars', action="store", dest='vcfTargets', type=str, help='vcf of loci to assess (posns only)', nargs='?', default=None)
parser.add_argument('-S','--snps', action="store_true", dest='snpsOnly', help='assess SNPs from vcf [default: indels only] ', default=False)
parser.add_argument('-I','--indels', action="store_true", dest='indelsOnly', help='assess indels from vcf [default: indels only] ', default=False)

parser.add_argument('-F','--flank', action="store", dest='flank', type=int, help='flank to parse out around refs (bp)', nargs='?', default=100)
parser.add_argument('-t','--tmp', action="store", dest='tmp', type=str, help='temp folder for calculations <tmp>', nargs='?', default='tmp')


args = parser.parse_args()

TMP="./"+args.tmp+"/"
if not path.exists(TMP):
    os.makedirs(TMP)

if args.vcfTargets is None:
    args.vcfTargets = args.vcfFull

#if neither vartype selected, chose indels only
if not args.snpsOnly and not args.indelsOnly:
    args.indelsOnly=True

print >>sys.stderr, args.vcfTargets
targetsFile = open(args.vcfTargets,'r')
targets=vcf.Reader(targetsFile)

name = path.basename(args.fasta_align).replace('.fasta','')

i=0
se = [] #set ends
sc = [] #set counts
#nb set counts used to find correct sequence in consensus fasta seqs

blockindex =  open(TMP+"var.intervals.idx","w")

#INDEXING DICTS: varname / interval with index / fileindex keys 
#                  - (prob a waste of memory, but fuck it)
blockI = dict()    #OVERALL INDEX
blockFP = dict()   #FILE POSITION INDEX
intervalsFP = dict()
intervalsI = dict()

#OUTPUT DICT
quals = dict()

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

chrlens =  _parse_chrlens_from_dict(args.fasta_cons)

for var in targets:
    #print var

    chrom = var.CHROM
    start = var.POS-args.flank
    if start <= 0: start = 1 
    end = var.POS + len(var.REF) +args.flank
    if end > chrlens[chrom]: end = chrlens[chrom]

    if args.indelsOnly and not var.is_indel:
        continue
    elif args.snpsOnly and not var.is_snp:
        continue

    i+=1
    telodist = chrlens[chrom]-var.POS
    if var.POS < telodist: telodist=var.POS
    
    if var.is_indel:
        altlen = 0
        for alt in var.ALT:
            if len(alt) > altlen: altlen = len(alt)
        varlen = len(var.REF) - altlen 

    elif var.is_snp:
        varlen=1
    vname = chrom+":"+str(var.POS)

    quals[vname,"L"] = varlen
    quals[vname,"TD"] = telodist


    #figure out which non-overlapping set to put it in:
    setfound = False
    s=0
    while not setfound:
#        print str(s)+":"+str(len(se))+" ",
        if s >= len(se):
#            print "+"
            se += [end]
            sc += [1]
            setfound = True
        elif se[s] < start-1:
#            print "S"
            se[s] = end
            sc[s] +=1
            setfound = True
        else:
            s+=1
#    print ""

#    print >>blocks,chrom+":"+str(start)+"-"+str(end)
    print >>blockindex,"\t".join([str(i),
                                 chrom+":"+str(start)+"-"+str(end),
                                 str(s),
                                 var.CHROM+":"+str(var.POS),
                                 var.REF+"/"+str(var.ALT)]
                                 )
    intervalsFP[(s,sc[s])] =  chrom+":"+str(start)+"-"+str(end)
    intervalsI[i] =           chrom+":"+str(start)+"-"+str(end)
    blockI[i] =               chrom+":"+str(var.POS)
    blockFP[(s,sc[s])] =      chrom+":"+str(var.POS)
#blocks.close()
blockindex.close()

#GET SIZE AND LEVENSHTEIN DISTANCE OF INTERVALS
fullVarsFile = open(args.vcfFull,'r')
fullVcf=vcf.Reader(fullVarsFile)
for i in intervalsI:
    # print >>sys.stderr,i, intervalsI[i]
    chrom, st, en = re.split("[:-]",intervalsI[i])
    vname = blockI[i]
#    print "\t".join([ intervalsI[i],"->",
#                      chrom, st, en])
    st = int(st); en = int(en)
    region=fullVcf.fetch(chrom, st, en)
    levdist = 0;
    dist=0;
    for v in region:
        if v.is_snp:
            levdist +=1
        elif v.is_deletion:
            levdist += len(v.REF)
            dist -= (len(v.REF)-1)
        elif v.is_indel:
            levdist += len(v.ALT[0])
            dist += (len(v.ALT[0])-1)
        else:
            print >>sys.stderr, "unknown var type found: "+str(v)
            print >>sys.stderr, "                    in: "+intervals[I]
    quals[(vname,'IL')] = (en-st)+1+dist
    quals[(vname,'LD')] = levdist
    

#PRINT INTO SEPARATE INTERVAL FILES
for s in range(0,len(se)):
    print >>sys.stderr,"writing interval set "+str(s)
    blocks = open(TMP+"var."+str(s)+".intervals","w")
    #get intervals for this set
    ints = [i for (s2,i) in intervalsFP if s2 == s] 
    for i in ints:
        print >>blocks,intervalsFP[(s,i)]
    blocks.close()


    print >>sys.stderr,"calculating consensus seq "+str(s)
    altref_command = ['java','-jar','/humgen/gsa-hpprojects/GATK/bin/GenomeAnalysisTK-3.4-0-g7e26428/GenomeAnalysisTK.jar',
                  '-T', 'FastaAlternateReferenceMaker',
                  '-R', args.fasta_cons,
                  '-o', TMP+'var.'+str(s)+'.intervals.fasta',
#                  '-L', chrom+":"+str(st+1)+"-"+str(en),
                  '-L', TMP+'var.'+str(s)+'.intervals',
                  '-V', args.vcfFull]
    print >>sys.stderr, " ".join(altref_command)
    check_output(altref_command)



print >>sys.stderr, args.fasta_align

#CONCATENATE FASTAS & SWAP IDS
fasta2 = open(TMP+'var.intervals.IDs.fasta','w')
for s in range(0,len(se)):
    fasta1 = open(TMP+'var.'+str(s)+'.intervals.fasta','r')
    for line in fasta1:
        header = re.search('^>(\d+)', line)
        if header is not None:
            ID = int(header.group(1))
            print >>fasta2, '>'+blockFP[(s,ID)]
        else:
            print >>fasta2,line,


tags=['AS', #alignment score (no matched bases?)
      'NM', #edit distance
      'XS', #suboptimal align score
      'XF', #support from forward/reverse
      'XE'  #no supporting seeds
      ]

bwa_command = ['bwa','bwasw',
               args.fasta_align, TMP+'var.intervals.IDs.fasta']
print >>sys.stderr, " ".join(bwa_command)
samfilename = TMP+name+'.sam'
#    samfilename = TMP+"realigns.sam"
realigns = open(samfilename,"w")
call(bwa_command,stdout=realigns)
realigns.close()


    #doesn't need sorting, should be in read order
    #pysam.sort("-S","-n", TMP+"realigns.sam", TMP+"realigns.sam")
samfile = pysam.AlignmentFile(samfilename, "r")

for a in samfile:
    #    print a.qname,
    #    print a.is_secondary,
    #    print a.qual
    q = str(a.qname)
    if (q,'N') not in quals:
        quals[(q,'N')]=0
        for tag in tags:
            quals[(q,tag)]=0
    quals[(q,'N')]+=1  # count No of aligns
    for tag,val in a.tags:
        quals[(q,tag)] += val

#print results
blocks = [blockI[b] for b in blockI]
# names = [n for b,n,t in quals]
# names = list(set(names))

print '#'+args.fasta_cons
print '#'+args.vcfFull

#print "block",
#for name in names:
#    print "\t".join([""] + [i+"."+name for i in ["N","AS","NM","XS"]]),

print "\t".join(["i","block","L","LD","N","AS","RS","NM","XS"])
#print ""
for b in blockI:
    block = blockI[b]
    print str(b)+"\t"+block,
    #for name in names:
    if (block, 'N') not in quals:
        print "\t".join(map(str,["",
                                 quals[(block,'L')],
                                 quals[(block,'TD')],
                                 quals[(block,'IL')],
                                 quals[(block,'LD')],
                                 0,
                                 0,
                                 0,
                                 0,
                                 0
                                 ])),
        
    else:
        maxscore = float(quals[(block,'IL')])
        score = round(int(quals[(block,'AS')])/maxscore,3)
        print "\t".join(map(str,["",
                                 quals[(block,'L')],
                                 quals[(block,'TD')],
                                 quals[(block,'IL')],
                                 quals[(block,'LD')],
                                 quals[(block,'N')],
                                 quals[(block,'AS')],
                                 score,
                                 quals[(block,'NM')],
                                 quals[(block,'XS')]
                                 ])),
    print "" 

exit(0)

#REMOVE EXIT COMMAND TO WRITE VCF (PROB NOT THAT USEFUL YET)

#write out new VCF
#reget full vcf
fullVarsFile = open(args.vcfFull,'r')
fullVcf=vcf.Reader(fullVarsFile)

vcfout = args.vcfFull
vcfout = vcfout.replace('vcf.gz','vcf')
vcfout = vcfout.replace('.vcf','.SW.vcf')
vcfoutF = open(vcfout,'w')
vcfwriter=vcf.Writer(vcfoutF,fullVcf)

for rec in fullVcf:
    block = rec.CHROM+":"+str(rec.POS)
    if block in blocks:
        maxscore = float(quals[(block,'IL')])
        rec.INFO['LD'] = quals[(block,'LD')]
        if (block,'AS') in quals:            
            score = round(int(quals[(block,'AS')])/maxscore,3)
            rec.INFO['NM'] = quals[(block,'NM')]
            rec.INFO['RS'] = score
        else:
            rec.INFO['NM'] = 0
            rec.INFO['RS'] = 0
        
    vcfwriter.write_record(rec)
vcfoutF.close()

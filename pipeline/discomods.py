#!/usr/bin/python

import sys
import re
import numpy as np
import h5py as h5

def parse_chrs_from_dict(fasta):
    seqdictF = fasta.replace(".fasta",".dict")
    print >>sys.stderr, seqdictF
    seqdict = open(seqdictF,'r')
    locs = []
    for line in seqdict:
        snres = re.search('SN:(\S+)', line)
        if snres is not None:
            seqname = snres.group(1)
            seqlen = re.search('LN:(\d+)', line).group(1)
#            print >>sys.stderr, seqname, seqlen
            locs += [(seqname,1,int(seqlen))]
    return locs

def parse_samples_from_tab(samplesFile):
    samples = dict()
    for line in open(samplesFile,'r'):
        F = line.strip().split("\t")
        [seqid, lane, dataset] = F[0:3]
        bamfile = F[4]
        samples[(seqid,lane,dataset)] = bamfile
    return samples

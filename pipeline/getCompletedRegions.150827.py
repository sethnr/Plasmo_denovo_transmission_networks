#!/usr/bin/python

import sys
import argparse
from string import *
from subprocess import call, check_output
import os.path as path
from math import floor,ceil
import re 
import numpy as np
import h5py as h5

from discomods import * 


parser = argparse.ArgumentParser(description='get completed regions from done_jobs file')


parser.add_argument('-s','--samples', action="store", dest='samples', type=str, help='samples file', nargs='?', default=None)
parser.add_argument('-f','--fasts', action="store", dest='refFasta', type=str, help='ref dict file', nargs='?', default=None)

parser.add_argument('-F','--File', action="store", dest='done', type=str, help='sample region file - eventually get from mongodb', default=None)
args = parser.parse_args()


chroms = parse_chrs_from_dict(args.refFasta)
samples = parse_samples_from_tab(args.samples)
samples = samples.keys()
samples.sort()


longest = max([e for c,s,e in chroms])
karyo = len(chroms)
sno = len(samples)

#make array of samples x chroms x chromlength
done = np.zeros((sno,karyo,longest+1),dtype=np.int)

ctoi = dict()
#set unreachable regions (before or after chrom) to -1
#keep 0 index as -1 so we can just use 1-based chrom pos
for i in range(0,karyo):
    c,s,e = chroms[i]
    ctoi[c]=i
    done[:,i,0]=-1
    done[:,i,e+1:]=-1

stoi = dict()
sample_names=[s+"_"+str(l) for s,l,d in samples] #ordered list for hpy5 file
for i in range(0,sno):
    s,l,d = samples[i]
    stoi[s+"_"+l]=i
    
#read sample_region file
doneRegions = open(args.done,'r')
doneSR=list()
for line in doneRegions:
#    print line;
    s,l,r = line.rstrip().split()
    s = s+'_'+str(l)
#    print r.split('\W')
    doneSR += [(s,r)]

for s,r in doneSR:
    c,st,en = re.split('\W',r)
    ci=ctoi[c]
    si=stoi[s]
    done[si,ci,int(st):(int(en)+1)]=1

#SAVE TO H5 FILE
f = h5.File("completed.hdf5", "w")
f.create_dataset("completed", data=done)
print sample_names
f.create_dataset("samples", data=sample_names)
f.create_dataset("chromosomes", data=[c for c,s,e in chroms])

print f
f.close()

## for ci in range(0,karyo):
##     c,s,e = chroms[ci]
##     doneCount = np.sum(np.equal(done[:,ci,:],1))
##     total = np.sum(np.greater(done[:,ci,:],-1))
##     print "\t".join(map(str,[c,doneCount,total,round(total/doneCount,5)]))

print "MB remaining:"
print "CHR",
for si in range(0,sno):
    s,l,d = samples[si]
    s = s+"_"+l
    print "\t"+s,
print ''
for ci in range(0,karyo):
    c,s,e = chroms[ci]
    print c,
    for si in range(0,sno):
        doneCount = np.sum(np.equal(done[si,ci,:],1))
        total = np.sum(np.greater(done[si,ci,:],-1))
#        print "\t"+str(round(total/doneCount,5)),
        remain = total-doneCount
        if remain > 1e5:
            print "\t"+str(round((total-doneCount)/1e6,2))+"m",
        elif remain > 1e2:
            print "\t"+str(round((total-doneCount)/1e3,2))+"k",
        else:
            print "\t"+str(round((total-doneCount),3)),
    print ''


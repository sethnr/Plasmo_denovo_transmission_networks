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


parser = argparse.ArgumentParser(description='get completed regions from done_jobs file')


parser.add_argument('-b','--blocksize', action="store", dest='blocksize', type=int, help='size of visualisation blocks', nargs='?', default=10000)
parser.add_argument('-H','--H5', action="store", dest='H5', type=str, help='h5 file for completed jobs', default="./completed.hdf5")
parser.add_argument('--bysample', action="store_true", dest='bySample', help='print completed for each sample in block?', default=False)
parser.add_argument('--include_inaccessible','-A', action="store_false", dest='accessible', help='include inaccessible regions', default=True)
args = parser.parse_args()


#READ H5 FILE
f = h5.File(args.H5, "r")
done = f['/completed/']
doable = f['/callable/']
samples = f['/samples/']
chroms = f['/chromosomes/']

#sno = len(samples)
#karyo=len(chroms)
sno,karyo,maxpos=done.shape
print done.shape
print sno
print karyo
print maxpos
print args.blocksize

if args.bySample:
    print >>sys.stderr, "CHR",
    for si in range(0,sno):
        s = samples[si]
        print >>sys.stderr,"\t"+s,
    print >>sys.stderr,''


for ci in range(0,karyo):
    c = chroms[ci]
    pi=-1
    print c+"\t",
    for si in range(0,sno):
        blocklen = np.sum(np.greater(done[si,ci,:],-1))
        
        notDone = 0;
        doablelen = np.sum(doable[ci,:])

        doneCount = np.sum(np.equal(done[si,ci,:],1))
        if blocklen>0:
            print "\t"+str(round(doneCount/float(blocklen),3)),
    print ''

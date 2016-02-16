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

blocks=0
accessible=0
inaccessible=0
incompleteBlocks=0
incompleteAccessible=0

print >>sys.stdout,"\t".join(["#chrom","start","end"]),
if args.bySample:
    for si in range(0,sno):
        s = samples[si]
        print >>sys.stderr,"\t"+s,
    print >>sys.stderr,''
if args.accessible:
    print >>sys.stdout,"\t".join(["accessable","acess_not_done"])
else:
    print >>sys.stdout,"\t".join(["accessable","all_not_done"])


for ci in range(0,karyo):
    c = chroms[ci]
    pi=-1
    print >>sys.stderr, "="*10+c+"="*10
    for ps in range(1,maxpos,args.blocksize):
        pi+=1
        pe=ps+args.blocksize-1
#        if pie > e: pie=e
        bigtotal = np.sum(np.greater(done[:,ci,ps:pe],-1))

        if bigtotal>0:
            print c+"\t",
            print "\t".join(map(str,[ps,pe])),
            blocks +=1
        

        notDone = 0;
        blocklen = np.sum(np.greater(done[1,ci,ps:pe+1],-1))
        doablelen = np.sum(doable[ci,ps:pe+1])

        if args.accessible:
            total = doablelen
        for si in range(0,sno):
            doneCount = np.sum(np.equal(done[si,ci,ps:pe+1],1))
            if not args.accessible:
                total = np.sum(np.greater(done[si,ci,ps:pe+1],-1))
            if total==0:
                break
#        print "\t"+str(round(total/doneCount,5)),
            remain = total-doneCount
            if args.bySample:
                if remain > 1e5:
                    print "\t"+str(round((total-doneCount)/1e6,2))+"m",
                elif remain > 1e2:
                    print "\t"+str(round((total-doneCount)/1e3,2))+"k",
                else:
                    print "\t"+str(round((total-doneCount),3)),
            if remain > 0:
                notDone += 1
        if bigtotal>0:
            print "\t"+str(round(doablelen/float(blocklen),3))+"\t"+str(notDone)
        if notDone > 0:
            incompleteBlocks += 1
            if round(doablelen/float(blocklen),2)==1:
                incompleteAccessible+=1
        if round(doablelen/float(blocklen),2)==1:
            accessible+=1
        elif blocklen > 0:
            inaccessible +=1
            
print >>sys.stderr,"\t".join(map(str,["blocks","incomplete","access","inaccessible","incomAcces"]))
print >>sys.stderr,"\t".join(map(str,[blocks,incompleteBlocks,accessible,inaccessible,incompleteAccessible]))

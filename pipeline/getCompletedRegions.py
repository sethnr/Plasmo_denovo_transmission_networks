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
parser.add_argument('-f','--fasta', action="store", dest='refFasta', type=str, help='ref dict file', nargs='?', default=None)
parser.add_argument('-N','--new', action="store_true", dest='newFile', help='setup new H5 file', default=False)
parser.add_argument('-H','--H5', action="store", dest='H5', type=str, help='H5 file file', nargs='?', default="./completed.hdf5")

parser.add_argument('-F','--File', action="store", dest='done', type=str, help='sample region file - eventually get from mongodb', default=None)
args = parser.parse_args()


samples = parse_samples_from_tab(args.samples)
samples = samples.keys()
samples.sort()
sno = len(samples)

sample_names = [s+"_"+l for s,l,d in samples]

chroms = parse_chrs_from_dict(args.refFasta)
longest = max([e for c,s,e in chroms])
karyo = len(chroms)



if not path.exists(args.H5):
    args.newFile=True

f = h5.File(args.H5, "a")
if args.newFile:
#make array of samples x chroms x chromlength
    print >>sys.stderr, "building empty HDF5 matrix"
#    callmat = np.ones((karyo,longest+1),dtype=np.int)

    #set up null array (-1 areas outside of chr length
#        print np.equal(done[0,i,:],0)
    if 'completed' in f:
        print >>sys.stderr, "dataset [/completed/] found; overwriting"
        del f['/completed/']
    done = np.zeros((sno,karyo,longest+1),dtype=np.int)
    f.create_dataset("completed", data=done)

    if 'callable' in f:
        print >>sys.stderr, "dataset [/callable/] found; overwriting"
        del f['/callable/']
    callmat = np.zeros((karyo,longest+1),dtype=np.int)
    f.create_dataset("callable", data=callmat)    

    if 'chromosomes' in f:
        print >>sys.stderr, "dataset [/chromosomes/] found; overwriting"
        del f['/chromosomes/']
    f.create_dataset("chromosomes", data=[c for c,s,e in chroms])

    if 'samples' in f:
        print >>sys.stderr, "dataset [/samples/] found; overwriting"
        del f['/samples/']
    f.create_dataset("samples", data=sample_names)

    done = f['/completed/']
    callmat = f['/callable/']
    samples = f['/samples/']

    for i in range(0,karyo):
        c,s,e = chroms[i]
        print >>sys.stderr, "initialising empty HDF5 matrix "+c
        done[:,i,0]=-1
        done[:,i,e+1:]=-1
#        callmat[i,0]=0
#        callmat[i,e+1:]=0

else:
    if 'completed' not in f:
        done = np.zeros((sno,karyo,longest+1),dtype=np.int)
        f.create_dataset("completed", data=done)
    if 'callable' not in f:
        callmat = np.zeros((karyo,longest+1),dtype=np.int)
        f.create_dataset("callable", data=callmat)    
    if 'chromosomes' not in f:
        f.create_dataset("chromosomes", data=[c for c,s,e in chroms])
    if 'samples' not in f:
        f.create_dataset("samples", data=sample_names)

    done = f['/completed/']
    callmat = f['/callable/']
    samples = f['/samples/']



#get c/s lookup tables
ctoi = dict()
for i in range(0,karyo):
    c,s,e = chroms[i]
    ctoi[c]=i

stoi = dict()
#sample_names=[s+"_"+str(l) for s,l,d in samples] #ordered list for hpy5 file
for i in range(0,sno):
    s = samples[i]
    stoi[s]=i
    
    
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


## for ci in range(0,karyo):
##     c,s,e = chroms[ci]
##     doneCount = np.sum(np.equal(done[:,ci,:],1))
##     total = np.sum(np.greater(done[:,ci,:],-1))
##     print "\t".join(map(str,[c,doneCount,total,round(total/doneCount,5)]))

print "MB remaining:"
print "CHR",
for si in range(0,sno):
    s = samples[si]
    print "\t"+s,
print ''
for ci in range(0,karyo):
    c,s,e = chroms[ci]
    print c,
    for si in range(0,sno):
        doneCount = np.sum(np.equal(done[si,ci,:],1))
        callablebase = np.equal(done[ci,:],1)
        callabledone = np.logical_and(np.equal(done[ci,:],1),np.equal(done[si,ci,:],1))
        total = np.sum(np.greater(done[si,ci,:],-1))
#        print "\t"+str(round(total/doneCount,5)),
#        remain = total-doneCount
        remain=np.sum(callablebase)-np.sum(callabledone)
        if remain > 1e5:
            print "\t"+str(round((total-doneCount)/1e6,2))+"m",
        elif remain > 1e2:
            print "\t"+str(round((total-doneCount)/1e3,2))+"k",
        else:
            print "\t"+str(round((total-doneCount),3)),
    print ''


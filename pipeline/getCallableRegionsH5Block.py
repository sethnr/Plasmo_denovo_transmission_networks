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

parser.add_argument('-f','--fasta', action="store", dest='refFasta', type=str, help='ref dict file', nargs='?', default=None)
parser.add_argument('-H','--H5', action="store", dest='H5', type=str, help='H5 file file', nargs='?', default="./completed.hdf5")
parser.add_argument('-c','--callable', action="store", dest='callable', type=str, help='callable file', default=None)
parser.add_argument('-N','--new', action="store_true", dest='newFile', help='setup new H5 file', default=False)
parser.add_argument('-s','--samples', action="store", dest='samples', type=str, help='samples file', nargs='?', default=None)
parser.add_argument('-B','--blocks', action="store_true", dest='block', help='block file (chr, start, end, .* , Callable [True/False]', default=False)

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
    done = np.zeros((sno,karyo,longest+1),dtype=np.int)
#    callmat = np.ones((karyo,longest+1),dtype=np.int)
    callmat = np.zeros((karyo,longest+1),dtype=np.int)

    #set up null array (-1 areas outside of chr length
#        print np.equal(done[0,i,:],0)
    if 'completed' in f:
        print >>sys.stderr, "dataset [/completed/] found; overwriting"
        del f['/completed/']
    f.create_dataset("completed", data=done)

    if 'callable' in f:
        print >>sys.stderr, "dataset [/callable/] found; overwriting"
        del f['/callable/']
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
        callmat[i,0]=0
        callmat[i,e+1:]=0

else:
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
callableFile = open(args.callable,'r')
i=0

print >>sys.stderr, "parsing callable file"

if args.block:
    for line in callableFile:
        i+=1
        F =  line.rstrip().split()
        chrom,st,en = F[:3]
        iscallable=(F[-1]=="True")
        st = int(st)
        en = int(en)
        ci = ctoi[chrom]
        if iscallable:
            callmat[ci,st:en]=1
else:
    c_callables=dict()
    for c,s,e in chroms:
        c_callables[c]=list()

    for line in callableFile:
        i+=1
        if (i % 1e6)==0: print >>sys.stderr, str(i/1e6)+"m parsed"
        chrom,pos,call = line.rstrip().split()
        pos = int(pos)
        call = int(call)
        c_callables[chrom] += [call]
        #done[:,ci,pos:(pos+1)]=-1
#        if int(call) ==1:
#            callmat[ci,pos:(pos+1)]=1
    for c,s,e in chroms:
        ci=ctoi[chrom]
        callmat[ci,s:e]=c_callables[c]
print done
print f


print "callable:"
for ci in range(0,karyo):
    c,s,e = chroms[ci]
    print c,ci,e
    print callmat[ci,:]
    callables = np.sum(callmat[ci,:])
    remain = e-callables
    callpc = callables/float(e)
    print "\t".join(map(str,[c,e,callables,callpc]))

f.close()

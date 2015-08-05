#!/usr/bin/python

import sys
import argparse
import subprocess
import os
from math import floor, ceil
import re
from string import *


parser = argparse.ArgumentParser(description='call STRs ')

#INPUTS:
parser.add_argument('-d','--depth', action="store", dest='dp', type=str, help='tabixed depth file for region', nargs='?', default=None)
parser.add_argument('-s','--swfile', action="store", dest='sw', type=str, help='tabixed sw file for region', nargs='?', default=None)
parser.add_argument('-r','--region','--regions', action="store", dest='region', type=str, help='region to call exhaustively across', nargs='?', default=None)
parser.add_argument('-b','--blocksize', action="store", dest='blocksize', type=int, default=1000, help='calc depth/SW in blocks of N')

parser.add_argument('-B','--base', action="store_true", dest='basewise', help='print out each base t/f value')

args = parser.parse_args()

locs=[]
if args.region is not None:
    if os.path.exists(args.region):
        regions = open(args.region,'r')
        for r in regions:
#            print >>sys.stderr, re.split('\W',r.rstrip())
            (rch,rst,ren) = re.split('\W',r.rstrip())
            locs += [(rch,int(rst),int(ren))]
    else:
        (rch,rst,ren) = re.split('\W',args.region)
        locs += [(rch,int(rst),int(ren))]

bl = args.blocksize
for ch,st,en in locs:
    for i in range(st,en,bl):
        ist = i+1
        ien = i+bl
        region = ch+":"+str(ist)+"-"+str(ien)
        swfile = open('tmp_sw.txt','w')
        dpfile = open('tmp_dp.txt','w')
        swcommand = ['tabix',args.sw,region]
        subprocess.call(swcommand,stdout=swfile)
        dpcommand = ['tabix',args.dp,region]
        subprocess.call(dpcommand,stdout=dpfile)
        swfile.close()
        dpfile.close()

        swfile = open('tmp_sw.txt','r')
        dpfile = open('tmp_dp.txt','r')

        NMpass = True
        RSpass = True
        
        NMlim = 0
        RSlim = 1
        RSfails=0
        NMfails=0
        NV=0
        
        varFailLim = 0.0


        for s in swfile:
            NV +=1
#            print len(s)
#            print s.split()
            i,c,p,block,IL,p2,L,LD,N,AS,RS,NM,XS = s.split()
            if int(NM) > NMlim: NMfails +=1
            if float(RS) > RSlim: RSfails +=1
        if NV > 0:
            if NMfails/float(NV) > varFailLim: NMpass=False
            if RSfails/float(NV) > varFailLim: RSpass=False
        
        #median = 54
        #stdev = 19.88
        dmax = 74
        dmin = 34

        dn = 0
        df = 0
        dFailLim = 0.1
        for d in dpfile:
            c,pos,dp = d.split()
            dn +=1
            if int(dp) > dmax: df +=1
            if int(dp) < dmin: df +=1
        Dpass=True
        if df/float(dn) > dFailLim: Dpass=False

        multiPass = (NMpass & RSpass & Dpass)
        if args.basewise:
            for i in range(ist,ien):
                print >>sys.stdout, "\t".join([ch,str(i),str(int(multiPass))])
        else:
            print >>sys.stdout, "\t".join(map(str,[ch,ist,ien,NV,NMfails,NMpass,RSfails,RSpass,df,dn,Dpass,multiPass]))

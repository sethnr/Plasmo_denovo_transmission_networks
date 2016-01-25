#!/bin/python

import sys
import argparse
from string import *
from os import path
import re 
parser = argparse.ArgumentParser(description='parse out statistics from VCF ')

parser.add_argument('-f','--freec', action="store", dest='freec', type=str, help='freec file1', nargs='?', default=None)
parser.add_argument('-pD','--pindelD', action="store", dest='pindelD', type=str, help='pindel deletions file', nargs='?', default=None)
parser.add_argument('-pTD','--pindelTD', action="store", dest='pindelTD', type=str, help='pindel tandem duplications file', nargs='?', default=None)
parser.add_argument('--mean', action="store_true", dest='mean', help='print mean distance & copy number',)
parser.add_argument('--min', action="store_true", dest='min', help='print min overlap & mean copy number',)
parser.add_argument('--max', action="store_true", dest='max', help='print combined distance & mean copy number',)

args = parser.parse_args()
freec = open(args.freec,'r')
pindelD = open(args.pindelD,'r')
pindelTD = open(args.pindelTD,'r')

#chrD,stD,enD,cnD,cnvD = None,None,None,None,None
#chrTD,stTD,enTD,cnTD,cnvTD = None,None,None,None,None

def readFREEC(cnvfile):
    line = cnvfile.readline()
#    print >>sys.stderr, cnvfile,"LINE:"+line+":",line is None,len(line),len(line)==0
    #if len(line)==0: line = None
#    if line is not None:
    if len(line)==0:
        (chr,st,en,cn,cnv) = (None,None,None,None,None)
    else:
        F = line.split()
        (chr,st,en,cn,cnv) = F[0:5]
        st = int(st); en=int(en); cn=int(cn)
    return (chr,st,en,cn,cnv)

def readPINDELS(fileD,fileTD):
    global chrD,stD,enD,cnD,cnvD
    global chrTD,stTD,enTD,cnTD,cnvTD
    #if last compared was TD, get new TD
    if ((chrP,stP,enP,cnP,cnvP) == (chrTD,stTD,enTD,cnTD,cnvTD)):
        (chrTD,stTD,enTD,cnTD,cnvTD) = readPINDEL(fileTD)
    #if last compared was D, get new D
    elif ((chrP,stP,enP,cnP,cnvP) == (chrD,stD,enD,cnD,cnvD)):
        (chrD,stD,enD,cnD,cnvD) = readPINDEL(fileD)
    #if one file has ended, return the other
    if chrD is None:
        return (chrTD,stTD,enTD,cnTD,cnvTD)
    elif chrTD is None:
        return (chrD,stD,enD,cnD,cnvD)
    else:
        return min((chrD,stD,enD,cnD,cnvD),(chrTD,stTD,enTD,cnTD,cnvTD))
    
def readPINDEL(pinfile):
#    print >>sys.stderr, cnvfile,"LINE:"+line+":",line is None,len(line),len(line)==0
    #if len(line)==0: line = None
#    if line is not None:
    line="INIT"; cn=1;cnv="."
    header=False

    line = pinfile.readline()
    while line:
        isheader=re.match('#', line)
#        print line, len(line),isheader
        if isheader:
#            print "header found"
            header=True
            line = pinfile.readline()
#            print line
            F = line.split()
#            print "\t".join(F[0:15])
            if F[1]=="D": cnv="loss"; cn=0
            elif F[1]=="TD": cnv="gain"; cn=2
#            (chr,st,en) = (F[7],F[9],F[10])
            (chr,st,en) = (F[7],F[12],F[13])   #use largest pindel range
            st = int(st); en=int(en); 
#            print "\t".join(map(str,[chr,st,en,(en-st),cnv,cn]))
            return (chr,st,en,cn,cnv)
        else:
            pass
#            print "no header"
        line = pinfile.readline()
    return (None,None,None,None,None)


#line1 = file1.readline()
#line2 = file2.readline()
#(chr1,st1,en1,cn1,cnv1) = line1.split()
#(chr2,st2,en2,cn2,cnv2) = line2.split()
#get tandem dupications and deletions, put earlies into pindel file 
(chrD,stD,enD,cnD,cnvD) = readPINDEL(pindelD)
(chrTD,stTD,enTD,cnTD,cnvTD) = readPINDEL(pindelTD)
(chrP,stP,enP,cnP,cnvP) = min((chrD,stD,enD,cnD,cnvD),(chrTD,stTD,enTD,cnTD,cnvTD))

(chrF,stF,enF,cnF,cnvF) = readFREEC(freec)

lastP=(None,None,None,None,None)
lastF=(None,None,None,None,None)

while (chrF is not None and chrP is not None):

    if 1==2:
        if (chrP,stP,enP,cnP,cnvP) == lastP:
            print "\t".join(map(str,[chrF,stF,enF,"   .   ","   .   ","   .   ",cnF,cnP,cnvF,cnvP]))
        elif (chrF,stF,enF,cnF,cnvF) == lastF:
            print "\t".join(map(str,["   .   ","   .   ","   .   ",chrP,stP,enP,cnF,cnP,cnvF,cnvP]))
        else:
            print "\t".join(map(str,[chrF,stF,enF,chrP,stP,enP,cnF,cnP,cnvF,cnvP]))

    if chrF == chrP and (enF > stP and stF < enP) and (cnvF==cnvP):
        pcOl = round((int(enF)-int(stF))/(float(enP)-int(stP)),2)        
#        print pcOl,
        minst = min(stF,stP)
        maxst = max(stF,stP)
        minen = min(enF,enP)
        maxen = max(enF,enP)
        meanst=(stF+stP)/2
        meanen=(enF+enP)/2
        meancn=(cnF+cnP)/2
        
        if args.mean:
            print "\t".join(map(str,[chrF,meanst,meanen,meancn,cnvF,pcOl]))
        elif args.max:
            print "\t".join(map(str,[chrF,minst,maxen,meancn,cnvF,pcOl]))
        elif args.min: 
            print "\t".join(map(str,[chrF,maxst,minen,meancn,cnvF,pcOl]))
        else:
            print "\t".join(map(str,[chrF,maxst,minen,minst,maxen,cnF,cnP,cnvF,pcOl]))
            print "\t".join([chrF,str(stF)+"-"+str(enF),str(enF-stF),str(stP)+"-"+str(enP),str(enP-stP),str(pcOl)])
    else:
#        print ".",
        pass

    lastP=(chrP,stP,enP,cnP,cnvP)
    lastF=(chrF,stF,enF,cnF,cnvF)
    
    if chrF > chrP:
#        print "F>P",
        (chrP,stP,enP,cnP,cnvP) = readPINDELS(pindelD,pindelTD)
    elif chrP > chrF:
#        print "P>F",
        (chrF,stF,enF,cnF,cnvF) = readFREEC(freec)
    elif chrF == chrP and (stF > stP or (stF==stP and enF>enP)):
#        print "F>P",
        (chrP,stP,enP,cnP,cnvP) = readPINDELS(pindelD,pindelTD)
    elif chrP == chrF and (stP > stF or (stF==stP and enP>enF)):
#        print "P>F",
        (chrF,stF,enF,cnF,cnvF) = readFREEC(freec)
    elif chrF == chrP and stF == stP and enF==enP:
#        print "P=F",
        (chrF,stF,enF,cnF,cnvF) = readFREEC(freec)
        (chrP,stP,enP,cnP,cnvP) = readPINDELS(pindelD,pindelTD)
#    print ""

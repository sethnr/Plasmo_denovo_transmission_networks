#!/usr/bin/python

#wrapper for discovar process
import vcf
import sys
import argparse
import subprocess
import re
from time import sleep
# import pythongrid as grid

from string import *

#SET VARS
RUNDISCO = '$DISCO1/runDiscovarVar.sh'



parser = argparse.ArgumentParser(description='call STRs ')

#INPUTS:
parser.add_argument('-t','--STR', action="store", dest='strLocs', type=str, help='bed file for STR locations', nargs='?', default=None)
parser.add_argument('-n','--SNP', action="store", dest='snpLocs', type=str, help='bed file for SNP locations', nargs='?', default=None)
parser.add_argument('-s','--samples', action="store", dest='samples', type=str, help='tab-delim file for sample names [sample_id, lane#, set]', nargs='?', default=None)
parser.add_argument('-d','--dir', action="store", dest='datadir', type=str, help='data directory for bam files', nargs='?', default=None)

#VARS
parser.add_argument('-z','--discosize', action="store", dest='discoSize', type=int, help='default size for expanded region for STR calling', nargs='?', default=5000)
parser.add_argument('-i','--defaultSize', action="store", dest='defSize', type=int, help='default size for STRs', nargs='?', default=0)
parser.add_argument('--dryrun', action="store_true", dest='dryrun', help='run through without bsubbing anything', default=False)


#DISCO INPUTS
parser.add_argument('-f','--ref', action="store", dest='refFasta', type=str, help='fasta file for reference seq', nargs='?', default=None)
parser.add_argument('-r','--region', action="store", dest='region', type=str, help='region to call exhaustively across', nargs='?', default=None)

args = parser.parse_args()


def _subjob(commandline, jobname, mem=2000, nodes=1, dependency=None):
    resource =  'select[mem>'+str(mem)+'] rusage[mem='+str(mem)+']  span[ptile='+str(nodes)+']'
    bsub_command = "-J "+jobname+" "+\
        "-e out/"+jobname+".e -o out/"+jobname+".o "+ \
        "-R "+resource+" -M "+str(mem)+" -n "+str(nodes) +" "+ \
        commandline

    print >>sys.stderr, "bsub",bsub_command
    bsub_return = """
Please specify a project.  You can set it with "bsub -P project ..."
or in the LSB_DEFAULTPROJECT environment variable.

By default, submitting under project "unspecified--broadfolk".
Job <1111> is submitted to queue <bhour>.
"""
    job_no = '-1'
        
    if not args.dryrun:
        bsub_return = subprocess.check_output("bsub", bsub_command, shell=True)   
        x = re.findall('Job <(\d+)> is submitted to queue <.*>',bsub_return)
        if len(x)>0:
            job_no = x[0]
    return job_no


def _waitdone(jobs):

    alldone = len(jobs)
    sleeptime = 0
    print "bjobs "+(" ".join(jobs))
    done = 0
    while done < alldone:
        done = 0
 
        sleep(sleeptime)
        bjobs_return = """JOBID   USER    STAT  QUEUE      FROM_HOST   EXEC_HOST   JOB_NAME   SUBMIT_TIME
    3789887 engreit DONE  hour       node1383    node1371    *uccessful May 11 17:12
    3789899 sredmon WAIT  bhour      tin         node1373    test       May 11 17:12"""
        if not args.dryrun:
            bjobs_return = subprocess.check_output("bjobs", " ".join(jobs), shell=True)
            for l in bjobs_return.split("\n"):
                F = l.strip().split()
                job_id, user, status = F[:3]
                if status == "DONE": done +=1
                    
        print >>sys.stderr, "DONE: ",done,"/",alldone
        if sleeptime < 60:
            sleeptime += 10
        if args.dryrun:
            done = alldone
            
samples = dict()
# parse tab file for samples / lanes / sets
for line in open(args.samples,'r'):
    F = line.strip().split("\t")
    [seqid, lane, dataset] = F[0:3]
    bamfile = "noBamCheck"
    # lsstring = args.datafile,"/",lane,"/",dataset,"/*bam"
    # bamfile = subprocess.check_output("ls", lsstring, shell=True)
    samples[(seqid,lane,dataset)] = bamfile
    print >>sys.stderr, seqid, lane, dataset, bamfile


flank = args.discoSize/2

SNPs = []
STRs = []
locs = []

# parse bed file SNPs
for line in open(args.snpLocs,'r'):
    F = line.strip().split("\t")
    [chrom, pos] = F[:2]
    pos = int(pos)
    st = (pos-flank); en = (pos+flank)

    SNPs += [(chrom,pos)]
    locs += [(chrom,st,en)]
    
# parse bed file STRs
for line in open(args.snpLocs,'r'):
    F = line.strip().split("\t")
    [chrom, start, end] = F[:3]
    start = int(start)
    end = int(end)
    st = (start-flank); en = (end+flank)

    STRs += [(chrom,start,end)]
    locs += [(chrom,st,en)]

#for (chrom1, st1, en1) in locs:
#    for (chrom2, st2, en2) in locs:

print locs

regions = ",".join([(chrom+":"+str(st)+"-"+str(en)) for (chrom,st,en) in locs])
print regions

# run discovar-var on bed file positions
# for all samples in tab-delim file
jobs = []
for (seqid, lane, dataset) in samples:
    command = " ".join([RUNDISCO, seqid, lane, regions])
    name = seqid+"_"+lane
    jobNo = _subjob(command, name )
    jobs += [jobNo]
    print seqid, lane, name, dataset, jobNo

_waitdone(jobs)

# make summary plots/files/graphs for samples in tab-delim file
# (subdivide into categories based on tab-delim file





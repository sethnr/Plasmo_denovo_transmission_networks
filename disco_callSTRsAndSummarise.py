#!/usr/bin/python

#wrapper for discovar process
import vcf
import sys
import argparse
import subprocess
import re
import os
from time import sleep
# import pythongrid as grid

from string import *

#SET VARS
RUNDISCO = '$DISCO1/runDiscovarVar.sh'
DATADIR = '/seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/'
RUNDIR = os.getcwd()

parser = argparse.ArgumentParser(description='call STRs ')

#INPUTS:
parser.add_argument('-t','--STR', action="store", dest='strLocs', type=str, help='bed file for STR locations', nargs='?', default=None)
parser.add_argument('-n','--SNP', action="store", dest='snpLocs', type=str, help='bed file for SNP locations', nargs='?', default=None)
parser.add_argument('-s','--samples', action="store", dest='samples', type=str, help='tab-delim file for sample names [sample_id, lane#, set]', nargs='?', default=None)
parser.add_argument('-d','--dir', action="store", dest='datadir', type=str, help='data directory for bam files', nargs='?', default=None)
parser.add_argument('--rundir', action="store", dest='rundir', type=str, help='directory for output (default = $PWD)', nargs='?', default=None)

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
    bsub_command = ["bsub",
                    "-J", jobname,
                    "-e", RUNDIR+"/out/"+jobname+".e",
                    "-o", RUNDIR+"/out/"+jobname+".o",
#                    "-e", jobname+".e",
#                    "-o", jobname+".o",
                    "-R", resource,
                    "-M", str(mem),
                    "-n", str(nodes),
                    commandline]

    print >>sys.stderr, "bsub",bsub_command
    bsub_return = """
Please specify a project.  You can set it with "bsub -P project ..."
or in the LSB_DEFAULTPROJECT environment variable.

By default, submitting under project "unspecified--broadfolk".
Job <1111> is submitted to queue <bhour>.
"""
    job_no = '-1'
    
    if not args.dryrun:
        bsub_return = subprocess.check_output(bsub_command)
        print >>sys.stderr, bsub_return
        
        x = re.findall('Job <(\d+)> is submitted to queue <.*>',bsub_return)
        if len(x)>0:
            job_no = x[0]
    return job_no


def _waitdone(jobs):

    alldone = len(jobs)
    sleeptime = 0
    print "bjobs "+(" ".join(jobs))
    done = 0
    bjobs_command = ["bjobs"] + jobs
    print >>sys.stderr, bjobs_command
    print >>sys.stderr, "FAIL\tPEND\tRUN\tDONE\tTOTAL"

    while done < alldone:
        done = 0; run = 0; pend = 0; fail = 0
        sleep(sleeptime)
        bjobs_return = """JOBID   USER    STAT  QUEUE      FROM_HOST   EXEC_HOST   JOB_NAME   SUBMIT_TIME
    3789887 engreit DONE  hour       node1383    node1371    *uccessful May 11 17:12
    3789899 sredmon WAIT  bhour      tin         node1373    test       May 11 17:12"""
        if not args.dryrun:
            bjobs_return = subprocess.check_output(bjobs_command)
        #    print >>sys.stderr, bjobs_return
            for l in bjobs_return.split("\n"):
                F = l.strip().split()
                if len(F) >=3:
                    job_id, user, status = F[:3]
                    if status == "DONE": done +=1
                    elif status == "PEND": pend +=1
                    elif status == "RUN": run += 1
                    elif status == "DIED": fail += 1
        sys.stderr.write(str(fail)+"\t"+str(pend)+"\t"+str(run)+"\t"+str(done)+"\t"+str(alldone)+"\r")
        sys.stderr.flush();
        
        if sleeptime < 60:
            sleeptime += 10
        if args.dryrun:
            done = alldone
    print >>sys.stderr, ""
            


#variables / defaults / etc
flank = args.discoSize/2
if args.datadir is not None: DATADIR = args.datadir
if args.rundir is not None:  RUNDIR = args.rundir



######
#PARSE FILES FOR INPUTS
######

# parse tab file for samples / lanes / sets
samples = dict()
for line in open(args.samples,'r'):
    F = line.strip().split("\t")
    [seqid, lane, dataset] = F[0:3]
    bamfile = "noBamCheck"
    # lsstring = args.datafile,"/",lane,"/",dataset,"/*bam"
    # bamfile = subprocess.check_output("ls", lsstring, shell=True)
    samples[(seqid,lane,dataset)] = bamfile
    print >>sys.stderr, seqid, lane, dataset, bamfile


SNPs = []
STRs = []
locs = []
# parse bed file SNPs
if args.snpLocs is not None:
    for line in open(args.snpLocs,'r'):
        F = line.strip().split("\t")
        [chrom, pos] = F[:2]
        pos = int(pos)
        st = (pos-flank); en = (pos+flank)

        SNPs += [(chrom,pos)]
        locs += [(chrom,st,en)]

# parse bed file STRs
if args.strLocs is not None:
    for line in open(args.strLocs,'r'):
        F = line.strip().split("\t")
        [chrom, start, end] = F[:3]
        start = int(start)
        end = int(end)
        st = (start-flank); en = (end+flank)

        STRs += [(chrom,start,end)]
        locs += [(chrom,st,en)]

#for (chrom1, st1, en1) in locs:
#    for (chrom2, st2, en2) in locs:

#consolidate locations:
#to write

regions = ",".join([(chrom+":"+str(st)+"-"+str(en)) for (chrom,st,en) in locs])
print regions



########
# RUN DISCOVAR
########

if args.rundir is not None:
    os.chdir(RUNDIR)
if not os.path.exists("./out"):
    os.mkdir("./out")

# run discovar-var on bed file positions
# for all samples in tab-delim file
jobs = []


for (seqid, lane, dataset) in samples:
    command = " ".join([RUNDISCO, seqid, lane, DATADIR, regions])
    name = seqid+"_"+lane
    jobNo = _subjob(command, name )
    jobs += [jobNo]
    print seqid, lane, name, dataset, jobNo

_waitdone(jobs)

########
# RUN ANALYSES
########


# make summary plots/files/graphs for samples in tab-delim file
# (subdivide into categories based on tab-delim file





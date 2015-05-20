#!/usr/bin/python

#wrapper for discovar process
import vcf
import sys
import argparse
import subprocess
import re
import os
from math import floor, ceil
from time import sleep
# import pythongrid as grid

from string import *

#SET VARS
RUNDISCO = '$DISCO1/runDiscovarVar.sh'
DATADIR = '/seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/'
RUNDIR = os.getcwd()

MEMORYGAP=1


parser = argparse.ArgumentParser(description='call STRs ')

#INPUTS:
parser.add_argument('-t','--STR', action="store", dest='strLocs', type=str, help='bed file for STR locations', nargs='?', default=None)
parser.add_argument('-n','--SNP', action="store", dest='snpLocs', type=str, help='bed file for SNP locations', nargs='?', default=None)
parser.add_argument('-s','--samples', action="store", dest='samples', type=str, help='tab-delim file for sample names [sample_id, lane#, set]', nargs='?', default=None)
parser.add_argument('-d','--dir', action="store", dest='datadir', type=str, help='data directory for bam files', nargs='?', default=None)
parser.add_argument('--rundir', action="store", dest='rundir', type=str, help='directory for output (default = $PWD)', nargs='?', default=None)

#VARS
parser.add_argument('-z','--discosize', action="store", dest='discoSize', type=int, help='default size for expanded region for STR calling', nargs='?', default=2000)
parser.add_argument('-i','--defaultSize', action="store", dest='defSize', type=int, help='default size for STRs', nargs='?', default=0)
parser.add_argument('--dryrun', action="store_true", dest='dryrun', help='run through sample parsing and folder creation without doing anything', default=False)
parser.add_argument('--discosucks', action="store_true", dest='nobsub', help='run through without bsubbing any disco jobs', default=False)
parser.add_argument('--shallow', action="store_true", dest='nodepth', help='run through sample parsing and folder creation without doing anything', default=False)


#DISCO INPUTS
parser.add_argument('-f','--ref', action="store", dest='refFasta', type=str, help='fasta file for reference seq', nargs='?', default=None)
parser.add_argument('-r','--region', action="store", dest='region', type=str, help='region to call exhaustively across', nargs='?', default=None)
parser.add_argument('--mem', action="store", dest='mem', type=int, help='memory to use for discovar [2000]', nargs='?', default=2000)
parser.add_argument('--nodes', action="store", dest='nodes', type=int, help='nodes to use for discovar [1]', nargs='?', default=1)
parser.add_argument('-Q','--queue', action="store", dest='queue', type=str, help='queue for LSF jobs [bhour]', nargs='?', default="bhour")
parser.add_argument('--splitregions', action="store_false", dest='oneregion', help='submit each region as separate job [false]')


args = parser.parse_args()


def _subjob(commandline, jobname, mem=2000, nodes=1, queue="bhour", dependency=None):
    resource =  'select[mem>'+str(mem)+'] rusage[mem='+str(mem)+']  span[ptile='+str(nodes)+']'
    bsub_command = ["bsub",
                    "-J", jobname,
                    "-e", RUNDIR+"/out/"+jobname+".e",
                    "-o", RUNDIR+"/out/"+jobname+".o",
#                    "-e", jobname+".e",
#                    "-o", jobname+".o",
                    "-q", queue,
                    "-R", resource,
                    "-M", str(mem),
                    "-n", str(nodes),
                    commandline]

#    print >>sys.stderr, ' '.join(bsub_command)
#    bsub_return = """
#Please specify a project.  You can set it with "bsub -P project ..."
#or in the LSB_DEFAULTPROJECT environment variable.
#
#By default, submitting under project "unspecified--broadfolk".
#Job <1111> is submitted to queue <bhour>.
#"""
    job_no = '-1'
    if not args.nobsub:
        FNULL = open(os.devnull, 'w')
        bsub_return = subprocess.check_output(bsub_command,stderr=FNULL)
        FNULL.close()
        
        x = re.findall('Job <(\d+)> is submitted to queue <.*>',bsub_return)
        if len(x)>0:
            job_no = x[0]
    return (job_no, ' '.join(bsub_command))


def _waitdone(jobs, faillog=sys.stderr):
    if type(jobs) is dict:
        pass
        #jobsd = copy(jobs)
        #jobs = [job for job in jobs]
    
    failing=set()
    alldone = len(jobs)
    sleeptime = 0
    print "bjobs "+(" ".join(jobs))
    
#    bjobs_command = ["bjobs"] + jobs
#    print >>sys.stderr, bjobs_command
    print >>sys.stderr, "FAIL\tPEND\tRUN\tDONE\tTOTAL"

    finished = 0
    done = 0
    fail = 0    
    while finished < alldone:
        running = 0 
        pend = 0 
        sleep(sleeptime)
#        bjobs_return = """JOBID   USER    STAT  QUEUE      FROM_HOST   EXEC_HOST   JOB_NAME   SUBMIT_TIME
#    3789887 engreit DONE  hour       node1383    node1371    *uccessful May 11 17:12
#    3789899 sredmon WAIT  bhour      tin         node1373    test       May 11 17:12"""
        if not args.nobsub:
            bjobs_command = ["bjobs"] + jobs
            bjobs_return = subprocess.check_output(bjobs_command)
        #    print >>sys.stderr, bjobs_return
            for l in bjobs_return.split("\n"):
                F = l.strip().split()
                if len(F) >=3:
                    job_id, user, status = F[:3]
                    if status == "DONE": 
                        done +=1
                        jobs.remove(job_id)
                    elif status == "PEND": pend +=1
                    elif status == "RUN": running += 1
                    elif status == "EXIT": 
                        fail += 1
                        jobs.remove(job_id)
                        if job_id not in failing:
                            if type(jobs) is dict:
                                print >>faillog, job_id, jobs[job_id]
                            else:
                                print >>faillog, job_id
                            failing.add(job_id)
                            
                            
        sys.stderr.write("\033[K")
        sys.stderr.write(" "+str(fail)+"\t"+str(pend)+"\t"+str(running)+"\t"+str(done)+"\t"+str(alldone)+"\r")
        sys.stderr.flush();
        
        #add failed to total
        finished = done + fail
        if sleeptime < 60:
            sleeptime += 10
        if args.nobsub:
            finished = alldone
    print >>sys.stderr, ""
    print >>sys.stderr, "----------------"
    print >>sys.stderr, " "+str(fail)+"\t"+str(pend)+"\t"+str(running)+"\t"+str(done)+"\t"+str(alldone)



# list of regions (chrom, start, end) join those within <joinflank> of each other
def _merge_regions(locs, joinflank=0):
    i=1
    locs = sorted(locs)
    while i < len(locs):
        #print i
        (c1, s1, e1) = locs[i-1]
        (c2, s2, e2) = locs[i]
        if (c1 == c2) and (e1 + joinflank > s2):
        #    print "removing",locs[i-1],locs[i]
            locs.remove(locs[i])
            locs.remove(locs[i-1])
            locs += [(c1, s1, e2)]
            locs = sorted(locs)
        else: 
            i += 1
    return locs

#variables / defaults / etc
flank = args.discoSize/2
if args.datadir is not None: DATADIR = args.datadir
if args.rundir is not None:  RUNDIR = args.rundir
if args.dryrun:
    args.nobsub = True
    args.nodepth = True


######
#PARSE FILES FOR INPUTS
######
samples = dict()

SNPs = []
STRs = []
locs = []

if args.region is not None:
    (rch,rst,ren) = re.split('\W',args.region)
    locs += [(rch,int(rst),int(ren))]


# parse tab file for samples / lanes / sets
for line in open(args.samples,'r'):
    F = line.strip().split("\t")
    [seqid, lane, dataset] = F[0:3]
    bamfile = "noBamCheck"
    # lsstring = args.datafile,"/",lane,"/",dataset,"/*bam"
    # bamfile = subprocess.check_output("ls", lsstring, shell=True)
    samples[(seqid,lane,dataset)] = bamfile
    print >>sys.stderr, seqid, lane, dataset, bamfile

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

#consolidate locations into larger regions:
if args.oneregion:
    regions = [",".join([(chrom+":"+str(st)+"-"+str(en)) for (chrom,st,en) in locs])]
else:
    regions = [(chrom+":"+str(st)+"-"+str(en)) for (chrom,st,en) in locs]
# print sorted(locs)

print >>sys.stderr, "merging regions within "+str(flank)
print >>sys.stderr, str(len(locs))+" --> ",
locs = _merge_regions(locs,flank)
print >>sys.stderr, len(locs)
print >>sys.stderr, locs
########
# RUN DISCOVAR
########

if args.rundir is not None:
    os.chdir(RUNDIR)
if not os.path.exists("./out"):
    os.mkdir("./out")

# run discovar-var on bed file positions
# for all samples in tab-delim file
jobs = dict()
failjobs = [] 

os.system("use Discovar")

datasets = set([dataset for (seqid, lane, dataset) in samples])
print >>sys.stderr, datasets

sublog = open("submitted.log",'w')
killscript = open("killAll.sh",'w')

for thisdataset in datasets:
    samplesInDataset = [(seqid, lane, dataset) for (seqid, lane, dataset) in samples if dataset == thisdataset]
    #make and change into dataset subdir
    if not os.path.exists("./"+thisdataset):
        os.mkdir("./"+thisdataset)
    os.chdir("./"+thisdataset)
   
#    discovars = " NUM_THREADS="+str(args.nodes)
#    discovars = " MAX_MEMORY_GB="+str(ceil(args.mem/1000))
    
    for (seqid, lane, dataset) in samplesInDataset:
        for region in regions:
            command = " ".join([RUNDISCO,
                                "-n", str(args.nodes),
                                "-m", str(int(ceil(args.mem/1000))-MEMORYGAP),
                                "-d", DATADIR,
                                seqid, lane, region])
            if args.oneregion:
                name = seqid+"_"+lane
            else:
                name = seqid+"_"+lane+"_"+region
            
            (jobNo, subd_command) = _subjob(command, name, args.mem, args.nodes, args.queue)
            #jobs += [jobNo]
            jobs[jobNo] = (seqid+"_"+str(lane), region)
            print >>sublog, '#',jobNo
            print >>sublog, subd_command
            print >>killscript,"bkill ",jobNo
            print seqid, lane, name, dataset, jobNo

    os.chdir("../")
sublog.close()
killscript.close()
os.chmod('killAll.sh',000755)

faillog = open("failed.log",'w',0)
        
_waitdone([job for job in jobs], faillog)

if len(failjobs) > 0:
#    faillog = open("failed.log",'w')
    for job in failjobs:
        print >>faillog, '#',job
        print >>faillog, jobs[job]

faillog.close()


########
# RUN ANALYSES
########

# make summary plots/files/graphs for samples in tab-delim file
# (subdivide into categories based on tab-delim file

for dataset in datasets:
    print >>sys.stderr, "getting bam depths ",dataset
    if not args.nodepth:
        for loc in [chrom+":"+str(st)+"-"+str(en) for (chrom,st,en) in locs]:
            os.system('samtools depth -r '+loc+' '+dataset+'/*/*bam > '+dataset+'_'+loc+'.depth')
        os.system('cat '+dataset+'_*.depth > '+dataset+'.depth')
        for loc in [chrom+":"+str(st)+"-"+str(en) for (chrom,st,en) in locs]:
            os.system('rm '+dataset+'_'+loc+'.depth')

        
    # parse bed file SNPs
    if args.snpLocs is not None:
        pass
        if not args.dryrun:
            SNPlocs = ",".join([chrom+":"+str(pos)+"-"+str(pos) for (chrom, pos) in SNPs])
            print >>sys.stderr," ".join(['bash',
                                         '../mergeVcfSNPs.sh',
                                         dataset,
                                         SNPlocs,
                                         dataset+'/*/*filtered.vcf.gz'
                                         ])
            subprocess.check_call(['bash',
                                   '../mergeVcfSNPs.sh',
                                    dataset,
                                    SNPlocs,
                                    dataset+'/*/*filtered.vcf.gz'
                                    ])

    # parse bed file STRs
    if args.strLocs is not None:
        if not args.dryrun:
            STRlocs = ",".join([chrom+":"+str(st)+"-"+str(en) for (chrom, st, en) in STRs])

            print >>sys.stderr," ".join(['bash',
                                         '../mergeVcfIndels.sh',
                                         dataset,
                                         STRlocs,
                                         dataset+'/*/*filtered.vcf.gz'
                                         ])
            subprocess.check_call(['bash',
                                   '../mergeVcfIndels.sh',
                                    dataset,
                                    STRlocs,
                                    dataset+'/*/*filtered.vcf.gz'
                                    ])
            subprocess.check_call(['python','../getSTRlengthFromVCF.py',
                                    '-v', dataset+'_STRs.vcf',
                                    '-o', dataset+'_STRs.calls'
                                    ])
        else: pass



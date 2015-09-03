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
import h5py as h5 

#from sub_lsf import *
#from sub_sge import *

#SET VARS
PIPELINE='/seq/plasmodium/sredmond/pfdisco/pipeline'
RUNDISCO=PIPELINE+'/runDiscovarVar.sh'
#RUNJOBARR=PIPELINE+'/runJobInArray.sh'

DATADIR = '/seq/picard/H2MGCBCXX/C1-516_2015-04-30_2015-05-03/'
RUNDIR = os.getcwd()
REFERENCE = "/seq/plasmodium/sredmond/refs/PlasmoDB-24_Pfalciparum3D7_Genome.fasta"
MEMORYGAP=1


parser = argparse.ArgumentParser(description='call STRs ')

#INPUTS:
parser.add_argument('-t','--STR', action="store", dest='strLocs', type=str, help='bed file for STR locations', nargs='?', default=None)
parser.add_argument('-n','--SNP', action="store", dest='snpLocs', type=str, help='bed file for SNP locations', nargs='?', default=None)
parser.add_argument('-s','--samples', action="store", dest='samples', type=str, help='tab-delim file for sample names [sample_id, lane#, set]', nargs='?', default=None)
parser.add_argument('-d','--dir', action="store", dest='datadir', type=str, help='data directory for bam files', nargs='?', default=None)
parser.add_argument('--rundir', action="store", dest='rundir', type=str, help='directory for output (default = $PWD)', nargs='?', default=None)
#parser.add_argument('--tmp', action="store", dest='tmpdir', type=str, help='directory for output (default = $PWD)', nargs='?', default=None)

parser.add_argument('--genome', action="store_true", dest='wholeGenome', help='run exhaustively across all regions in genome', default=False)

parser.add_argument('--prefix', action="store", dest='prefix', help='prefix for submission log and killscripts (to allow concurrent runs in same folder', default="")

#VARS
parser.add_argument('-z','--discosize', action="store", dest='discoSize', type=int, help='default size for expanded region for STR calling', nargs='?', default=2000)
parser.add_argument('-i','--defaultSize', action="store", dest='defSize', type=int, help='default size for STRs', nargs='?', default=0)
parser.add_argument('--dryrun', action="store_true", dest='dryrun', help='run through sample parsing and folder creation without doing anything', default=False)
parser.add_argument('--discosucks', action="store_true", dest='nobsub', help='run through without bsubbing any disco jobs', default=False)
parser.add_argument('--shallow', action="store_true", dest='nodepth', help='run through sample parsing and folder creation without doing anything', default=False)
parser.add_argument('--interactive', action="store_true", dest='interactive', help='monitor jobs as being run (shell must remain open to finish analysis), default=[fire-and-forget]', default=False)

parser.add_argument('--donejobs', action="store", dest='donejobs', help='list of jobs that have already completed (remove from jobarray before submission)', default=None)
parser.add_argument('--doneH5', action="store", dest='doneH5', help='h5file of regions/samples that have already completed (if whole of job is covered, do not submit)', default=None)


#DISCO INPUTS
parser.add_argument('-f','--ref', action="store", dest='refFasta', type=str, help='fasta file for reference seq', nargs='?', default=None)
parser.add_argument('-r','--region','--regions', action="store", dest='region', type=str, help='region to call exhaustively across', nargs='?', default=None)
parser.add_argument('--mem', action="store", dest='mem', type=int, help='memory to use for discovar [2000]', nargs='?', default=2000)
parser.add_argument('--nodes', action="store", dest='nodes', type=int, help='nodes to use for discovar [1]', nargs='?', default=1)
parser.add_argument('--maxconc', action="store", dest='maxconc', type=int, help='maximum number of discovar jobs to run at once [20]', nargs='?', default=20)
parser.add_argument('-Q','--queue', action="store", dest='queue', type=str, help='queue for LSF jobs [bhour]', nargs='?', default="short")
parser.add_argument('--splitregions', action="store_false", dest='oneregion', help='submit each region as separate job [false]')
parser.add_argument('--splitsize', action="store", dest='splitsize', type=int, default=-1, help='split regions larger than N [-1]')
parser.add_argument('--splitflank', action="store", dest='splitflank', type=int, default=500, help='overlap regions by flank N bp [500]')
parser.add_argument('--minsplit', action="store", dest='minsplit', type=int, default=2000, help='do not split into chunk of less than N bp [2000]')

parser.add_argument('--bam_bylane', action="store_true", dest='bylane', default=True, help='parse bam location from lane [True]')
parser.add_argument('--bam_direct', action="store_false", dest='bylane', default=True, help='parse bam location from lane [True]')

parser.add_argument('-W','--jobtime', action="store", dest='jobtime', type=str, help='time for job', nargs='?', default=None)

parser.add_argument('--lsf', action="store_true", dest='lsf', default=False, help='sub with lsf rather than sge')

args = parser.parse_args()

if args.lsf:
    from sub_lsf import *
else: 
    from sub_sge import *





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

def _split_regions(locs, ssize=50000, sflank=500):
    locs = sorted(locs)
    #nb need to copy array, otherwise keeps iterating through new ranges
    for bigloc in locs[:]:
        (c1, s1, e1) = bigloc
        if e1 - s1 > ssize:
            #print "removing",bigloc
            if args.oneregion is True:
                args.oneregion = False
                print >> sys.stderr, "splitting large region, submitting separate jobs for each locus"
                
            locs.remove(bigloc)
            #print >>sys.stderr, range(s1, e1, ssize)
            for s2 in range(s1, e1, ssize):
                e2 = s2 + ssize
                if e2 > e1: e2 = e1
                if s2 > s1: s2 -= sflank
                if e2-s2 < args.minsplit:
                    s2=e2-args.minsplit
                locs += [(c1, s2, e2)]
                #print "adding", (c1, s2, e2)
                
            locs = sorted(locs)
    return locs

def _parse_chrs_from_dict(fasta):
    seqdictF = fasta.replace(".fasta",".dict")
    print >>sys.stderr, seqdictF
    seqdict = open(seqdictF,'r')
    locs = []
    for line in seqdict:
        snres = re.search('SN:(\S+)', line)
        if snres is not None:
            seqname = snres.group(1)
            seqlen = re.search('LN:(\d+)', line).group(1)
#            print >>sys.stderr, seqname, seqlen
            locs += [(seqname,1,int(seqlen))]
    return locs



#variables / defaults / etc
flank = args.discoSize/2
if args.donejobs is not None:
    args.donejobs = os.path.abspath(args.donejobs)
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


if args.wholeGenome and args.region is not None:
    print >>sys.stderr, "--genome and --region cannot both be set"
    exit(1)
elif args.wholeGenome:
    locs = _parse_chrs_from_dict(args.refFasta)
elif args.region is not None:
    if os.path.exists(args.region):
        regions = open(args.region,'r')
        for r in regions:
#            print >>sys.stderr, re.split('\W',r.rstrip())
            (rch,rst,ren) = re.split('\W',r.rstrip())
            locs += [(rch,int(rst),int(ren))]
    else:
        (rch,rst,ren) = re.split('\W',args.region)
        locs += [(rch,int(rst),int(ren))]

# parse tab file for samples / lanes / sets
for line in open(args.samples,'r'):
    F = line.strip().split("\t")
    [seqid, lane, dataset] = F[0:3]
    bamfile = F[4]
    samples[(seqid,lane,dataset)] = bamfile

# CHECK SHAPE OF SAMPLE FILE? DISTINGUISH BETWEEN sample/lane/dataset or name/bamfile
# LATER...
#    if len(F)==3:
#        [name, dataset, bamfile] = F[0:3]
#        samples[(name, dataset)] = bamfile
#    elif len(F)==4:
#        [seqid, lane, dataset, bamfile] = F[0:4]
#        samples[(seqid,lane,dataset,bamfile)] = bamfile
    
    # lsstring = args.datafile,"/",lane,"/",dataset,"/*bam"
    # bamfile = subprocess.check_output("ls", lsstring, shell=True)
    
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

#print >>sys.stderr, "merging regions within "+str(flank)
#print >>sys.stderr, str(len(locs))+" --> ",
locs = _merge_regions(locs,flank)
if args.splitsize > 0:
    locs = _split_regions(locs,args.splitsize,args.splitflank)
            
#consolidate locations into larger regions:
if args.oneregion:
    regions = [",".join([(chrom+":"+str(st)+"-"+str(en)) for (chrom,st,en) in locs])]
else:
    regions = [(chrom+":"+str(st)+"-"+str(en)) for (chrom,st,en) in locs]
# print sorted(locs)

if args.refFasta is not None:
    REFERENCE = args.refFasta

#print >>sys.stderr, locs
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

# os.system("use Discovar")

datasets = set([dataset for (seqid, lane, dataset) in samples])
print >>sys.stderr, datasets

print >>sys.stderr, "assessing "+str(len(locs))+" regions in "+str(len(samples))+" samples"
 
#commands = []
commands = dict()
#CONSTRUCT ALL COMMANDS:
for thisdataset in datasets:
    samplesInDataset = [(seqid, lane, dataset, samples[(seqid, lane, dataset)]) for (seqid, lane, dataset) in samples if dataset == thisdataset]

    for (seqid, lane, dataset, bamfile) in samplesInDataset:
        for region in regions:
            if args.bylane is True:
                command = " ".join([RUNDISCO,
                                 "-n", str(args.nodes),
                                 "-m", str(int(ceil(args.mem/1000))-MEMORYGAP),
                                 "-f", REFERENCE,
                                 "-d", DATADIR,
                                 seqid, lane, region])
            else:
                command = " ".join([RUNDISCO,
                                "-n", str(args.nodes),
                                "-m", str(int(ceil(args.mem/1000))-MEMORYGAP),
                                "-f", REFERENCE,
                                "-N", seqid,
                                "-B", bamfile,
                                    region])
            #commands += [command]
            commands[(seqid+'_'+lane,region)]=command
    


    jobname = thisdataset+"_disco"

    
    totalCommands = len(commands)

    if args.donejobs is not None:
        commands = commands.values()
        done = open(args.donejobs,'r')
        killlist = done.readlines()
        print >>sys.stderr, "killlist found: "+args.donejobs
        print >>sys.stderr, commands[1]
        print >>sys.stderr, killlist[1]
        print >>sys.stderr, "commands: "+str(len(commands))+" kills: "+str(len(killlist)),
        killlist = [k.rstrip() for k in killlist]
        commands = set(commands)-set(killlist)
        #commands = [c for c in commands if c not in killlist]
        print >>sys.stderr, " remain: "+str(len(commands))       

    elif args.doneH5 is not None:
        h5f= h5.File(args.doneH5, "r")

        h5samples = h5f['samples']
        stoi = dict()
        for i in range(0,len(h5samples)):
            stoi[h5samples[i]]=i

        h5chromosomes = h5f['chromosomes']
        ctoi = dict()
        for i in range(0,len(h5chromosomes)):
            ctoi[h5chromosomes[i]]=i

        callableJobs = dict()
        uncallable = list()
        alreadyDone=list()
        toDo=list()

        callmat = h5f['callable']
        done = h5f['completed']
        
        minCallablePC = 0.5
        print >>sys.stderr, "removing inacessible regions"
        for smp, region in commands:
            c,s,e = re.split('\W',region)
            s = int(s); e = int(e)
            iscallable = callmat[ctoi[c],s:e]
            jobdone = done[stoi[smp],ctoi[c],s:e]
            if sum(iscallable) <= (e-s)*minCallablePC:
                uncallable += [commands[(smp,region)]]
                pass
            else:
                callableJobs[(smp,region)] = commands[(smp,region)]
        commands = callableJobs

        print >>sys.stderr, "removing called jobs"
        for smp, region in commands:
            c,s,e = re.split('\W',region)
            s = int(s); e = int(e)
            jobdone = done[stoi[smp],ctoi[c],s:e]
            if sum(jobdone) == len(jobdone):
                alreadyDone += [commands[(smp,region)]]
                pass
            else:
                toDo += [commands[(smp,region)]]
        print >>sys.stderr, "commands: "+str(totalCommands), 
        print >>sys.stderr, "inaccess: "+str(len(uncallable)),
        print >>sys.stderr, "    done: "+str(len(alreadyDone)),
        print >>sys.stderr, "  remain: "+str(len(toDo)),

        commands = toDo
    else:
        commands = commands.values()

    #make and change into dataset subdir
    if not os.path.exists("./"+thisdataset):
        os.mkdir("./"+thisdataset)
    os.chdir("./"+thisdataset)
   
    
    (jobNo, subdcommand) = subarray(commands, jobname, args.mem, args.nodes, args.queue, args.jobtime, maxconc=args.maxconc, nobsub=args.nobsub)
    #return jobNo
    os.chdir("../")

    joblog(jobNo, subdcommand, args.prefix)
    
#    print >>sublog, '#',jobNo
#    print >>sublog, repr(subdcommand)
#    print >>killscript,"bkill ",jobNo
#    #print seqid, lane, name, dataset, jobNo
    jobs[thisdataset] = jobNo
    #os.chdir("../")
#    sublog.close()
#    killscript.close()
#os.chmod(killscriptname,000755)


#if running interactively, wait for jobs to finish and monitor
if args.interactive is True:
    faillog = open("failed.log",'w',0)
    _waitdone([jobs[dataset] for dataset in jobs], faillog, nobsub=args.nobsub)
    if len(failjobs) > 0:
        for job in failjobs:
            print >>faillog, '#',job
            print >>faillog, jobs[job]
    faillog.close()


########
# RUN ANALYSES
########

# make summary plots/files/graphs for samples in tab-delim file
# (subdivide into categories based on tab-delim file

#reopen log files
#sublog = open(sublogname,'a')
#killscript = open(killscriptname,'a')

for dataset in datasets:
    print >>sys.stderr, "merge/concat ",dataset
#    discoFinished="done("+jobs[dataset]+")"
    discoFinished=jobs[dataset]
    
    if args.interactive is True:
        # following code block for interactive running, for jobarray skip to next phase
        subprocess.check_call([
#                               'bash',
                               PIPELINE+'/mergeConcatVcfs.sh',
                               dataset])
    else:
        mergeCommand = " ".join([
#                                 'bash',
                                 PIPELINE+'/mergeConcatVcfs.sh',
                                 dataset])
        name = dataset+"_concat";
        (jobNo, subdcommand) = subjob(mergeCommand, name, args.mem, 1, args.queue, args.jobtime, dependency=discoFinished,nobsub=args.nobsub)
        ## print >>sublog, '#',jobNo
        ## print >>sublog, repr(subdcommand)
        ## print >>killscript,"bkill ",jobNo
        print >>sys.stderr,"MERGECOMMAND: "+mergeCommand
        print >>sys.stderr,"SUB'DCOMMAND: "+subdcommand
        print >>sys.stderr,"PWD: "+os.getcwd()
        
        joblog(jobNo, subdcommand, args.prefix) 




####
# remove all below? 
####

                
for dataset in datasets:
    print >>sys.stderr, "getting bam depths ",dataset
#        discoFinished="\"done("+jobs[dataset]+")\""
    discoFinished=jobs[dataset]

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
                                         PIPELINE+'/mergeVcfSNPs.sh',
                                         dataset,
                                         SNPlocs,
                                         dataset+'/*/*filtered.vcf.gz'
                                         ])
            if args.interactive is True:
                # following code block for interactive running, for jobarray skip to next phase
                subprocess.check_call(['bash',
                                   PIPELINE+'/mergeVcfSNPs.sh',
                                    dataset,
                                    SNPlocs,
                                    dataset+'/*/*filtered.vcf.gz'
                                    ])
            else:
                SNPcommand = " ".join(['bash',
                                   PIPELINE+'/mergeVcfSNPs.sh',
                                    dataset,
                                    SNPlocs,
                                    dataset+'/*/*filtered.vcf.gz'
                                    ])
                name = dataset+"_mergeSNPs";
                (jobNo, subdcommand) = subjob(SNPcommand, name, args.mem, args.nodes, args.queue, args.jobtime, dependency=discoFinished,nobsub=args.nobsub)
                ## print >>sublog, '#',jobNo
                ## print >>sublog, repr(subdcommand)
                ## print >>killscript,"bkill ",jobNo
                joblog(jobNo, subdcommand, args.prefix) 
            
    # Parse bed file STRs
    if args.strLocs is not None:
        if not args.dryrun:
            STRlocs = ",".join([chrom+":"+str(st)+"-"+str(en) for (chrom, st, en) in STRs])

            print >>sys.stderr," ".join([
#                                         'bash',
                                         PIPELINE+'/mergeVcfIndels.sh',
                                         dataset,
                                         STRlocs,
                                         dataset+'/*/*filtered.vcf.gz'
                                         ])
            if args.interactive is True:
                #write an interactive version soon!
                pass
            else:

                STRcommand1 = " ".join([
#                                        'bash',
                                        PIPELINE+'/mergeVcfIndels.sh',
                                        dataset,
                                        STRlocs,
                                        dataset+'/*/*filtered.vcf.gz'
                                        ])
            
                name = dataset+"_mergeSTRs";
             
                (jobNo1, subdcommand) = subjob(STRcommand1, name, args.mem, args.nodes, args.queue, args.jobtime, dependency=discoFinished)
                ## print >>sublog, '#',jobNo1
                ## print >>sublog, repr(subdcommand)
                ## print >>killscript,"bkill ",jobNo1
    
                subjob(jobNo1)
                STRcommand2 = " ".join(['python',PIPELINE+'/getSTRlengthFromVCF.py',
                                    '-v', dataset+'_STRs.vcf',
                                    '-o', dataset+'_STRs.calls'
                                    ])
                (jobNo, subdcommand) = subjob(STRcommand2, name, args.mem, args.nodes, args.queue, args.jobtime, dependency=jobNo2,nobsub=args.nobsub)
                ## print >>sublog, '#',jobNo2
                ## print >>sublog, repr(subdcommand)
                ## print >>killscript,"bkill ",jobNo2
                subjob(jobNo)
            
        else: pass



#!/usr/bin/python

import sys
import os 
import re
import subprocess

PIPELINE='/seq/plasmodium/sredmond/pfdisco/pipeline'
RUNJOBARR=PIPELINE+'/runJobInArray.sh'
RUNDIR = os.getcwd()

def doSomeStuff():
    print >>sys.stderr, "did some stuff"

def subarray(joblist, jobname, mem=2000, nodes=1, queue="bhour", jobtime=None, dependency=None, maxconc=-1):
    if dependency is not None:
        dependency="\"done("+dependency+")\""
    comlog = "sub_commands.txt"
    print >>sys.stderr, "OPENING COMLOG"
    comlogf = open(comlog,'w')
    
    for command in joblist:
        print >>sys.stderr, command
        print >>comlogf, command
    
    print >>sys.stderr, "CLOSING COMLOG"
    comlogf.close()
#    os.system("cp sub_commands.txt sub_commands.log")
#    os.system("wc -l sub_commands.*")
    arrCommand=RUNJOBARR+" "+comlog
    
    jobarrname = jobname+"["+str(1)+"-"+str(len(joblist))+"]"
    outname = jobname+".%I"
    if maxconc > 0: jobarrname += "%"+str(maxconc)
    (jobNo, subdcommand) = subjob(arrCommand, jobarrname, mem, nodes, queue, jobtime, dependency, outname)
    return (jobNo, subdcommand)
#    (jobNo, subdcommand) = subjob(commandline, jobname, mem=2000, nodes=1, queue="bhour", jobtime, dependency):


def subjob(commandline, jobname, mem=2000, nodes=1, queue="bhour", jobtime=None, dependency=None, outname=None, nobsub=False):
    if dependency is not None:
        dependency="\"done("+dependency+")\""
    resource =  'select[mem>'+str(mem)+'] rusage[mem='+str(mem)+']  span[ptile='+str(nodes)+']'
    if outname is None: outname=jobname
    bsub_command = ["bsub",
                    "-J", jobname,
                    "-e", RUNDIR+"/out/"+outname+".e",
                    "-o", RUNDIR+"/out/"+outname+".o",
#                    "-e", jobname+".e",
#                    "-o", jobname+".o",
                    "-q", queue,
                    "-R", resource,
                    "-M", str(mem),
                    "-n", str(nodes)]
    if dependency is not None: bsub_command += ['-w',dependency]
    if jobtime is not None: bsub_command += ['-W',jobtime]
    bsub_command += [commandline]
    
#    print >>sys.stderr, ' '.join(bsub_command)
#    bsub_return = """
#Please specify a project.  You can set it with "bsub -P project ..."
#or in the LSB_DEFAULTPROJECT environment variable.
#
#By default, submitting under project "unspecified--broadfolk".
#Job <1111> is submitted to queue <bhour>.
#"""
    job_no = '-1'
    if not nobsub:
        FNULL = open(os.devnull, 'w')
        bsub_return = subprocess.check_output(bsub_command,stderr=FNULL)
        FNULL.close()
        
        x = re.findall('Job <(\d+)> is submitted to queue <.*>',bsub_return)
        if len(x)>0:
            job_no = x[0]
    return (job_no, ' '.join(bsub_command))


def waitdone(jobs, faillog=sys.stderr, nobsub=False):
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
        if not nobsub:
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
        if nobsub:
            finished = alldone
    print >>sys.stderr, ""
    print >>sys.stderr, "----------------"
    print >>sys.stderr, " "+str(fail)+"\t"+str(pend)+"\t"+str(running)+"\t"+str(done)+"\t"+str(alldone)


def joblog(job_id, subdcommand, prefix=None, sublogname="submitted.log", killscriptname="killAll.sh"):
    if prefix is not None:
        sublogname = prefix+"."+sublogname
        killscriptname = prefix+"."+killscriptname

    if not os.path.exists(killscriptname):
        killscript = open(killscriptname,'a')
        os.chmod(killscriptname,000755)
    else:        
        killscript = open(killscriptname,'a')
    sublog = open(sublogname,'a')

    print >>sublog, '#',job_id
    print >>sublog, repr(subdcommand)
    print >>killscript,"bkill ",job_id
    sublog.close()
    killscript.close()
    

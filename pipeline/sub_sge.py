#!/usr/bin/python

import sys
import os 
import re
import subprocess
from math import ceil,floor

PIPELINE='/seq/plasmodium/sredmond/pfdisco/pipeline'
RUNJOBARR=PIPELINE+'/runJobInArray.sh'
RUNDIR = os.getcwd()

def subarray(joblist, jobname, mem=2000, nodes=1, queue="bhour", jobtime=None, dependency=None, maxconc=-1, nobsub=False):
    comlog = "sub_commands.txt"
    print >>sys.stderr, "OPENING COMLOG"
    print >>sys.stderr, "dryrun",nobsub
    comlogf = open(comlog,'w')
    
    for command in joblist:
        #print >>sys.stderr, command
        print >>comlogf, command
    
    print >>sys.stderr, "CLOSING COMLOG"
    comlogf.close()
#    os.system("cp sub_commands.txt sub_commands.log")
#    os.system("wc -l sub_commands.*")

    print >>sys.stderr, comlog,"\n",os.path.abspath(comlog)
    arrCommand=RUNJOBARR+" "+comlog
    
#    jobarrname = jobname+"["+str(1)+"-"+str(len(joblist))+"]"
    arrayjobs=len(joblist)
    outname = jobname
#    if maxconc > 0: jobarrname += "%"+str(maxconc)
    (jobNo, subdcommand) = subjob(arrCommand, jobname, mem, nodes, queue, jobtime, dependency, outname, arrayjobs, maxconc,nobsub)
    return (jobNo, subdcommand)
#    (jobNo, subdcommand) = subjob(commandline, jobname, mem=2000, nodes=1, queue="bhour", jobtime, dependency):


def subjob(commandline, jobname, mem=2000, nodes=1, queue="bhour", jobtime=None, dependency=None, outname=None, arrayjobs=-1, maxconc=-1, nobsub=False):
#    resource =  'select[mem>'+str(mem)+'] rusage[mem='+str(mem)+']  span[ptile='+str(nodes)+']'
    resource =  [
                 '-l','m_mem_free='+str(floor(mem/1000))+'G',
#                 '-l','mem='+str(floor(mem/1000))+'G',
        
                 '-l','m_core='+str(nodes)]                        
    if jobtime is not None: resource += ['-l','h_rt='+jobtime]

    if outname is None: outname=jobname
#    qsub_command = ["qsub",
#                    "-J", jobname,
#                    "-e", RUNDIR+"/out/"+outname+".e",
#                    "-o", RUNDIR+"/out/"+outname+".o",
#                    "-q", queue,
#                    "-R", resource,
#                    "-M", str(mem),
#                    "-n", str(nodes)]
    if re.match('\d',jobname[0]): jobname = 'j'+jobname
    qsub_command = ["qsub",'-terse',
                    "-N", jobname,
#                    "-e", RUNDIR+"/out/"+outname+".${SGE_TASK_ID}.e",
#                    "-o", RUNDIR+"/out/"+outname+".${SGE_TASK_ID}.o",
                    "-e", RUNDIR+"/out/",
                    "-o", RUNDIR+"/out/",
                    "-q", queue]
    qsub_command += resource
    if arrayjobs > 0: qsub_command +=  ['-t','1-'+str(arrayjobs)]
    if maxconc > 0: qsub_command += ['-tc',str(maxconc)]

#    -W depend=afterok:<Job-ID>
    
    if dependency is not None:
#        dependency="depend=afterok:"+dependency
#        qsub_command += ['-W',dependency]
        qsub_command += ['-hold_jid',dependency]
    qsub_command += commandline.split()
    
#    print >>sys.stderr, ' '.join(qsub_command)
#    qsub_return = """
#Please specify a project.  You can set it with "qsub -P project ..."
#or in the LSB_DEFAULTPROJECT environment variable.
#
#By default, submitting under project "unspecified--broadfolk".
#Job <1111> is submitted to queue <bhour>.
#"""

    print >>sys.stderr, " ".join(qsub_command)

    job_no = '-1'
    if not nobsub:
#        FNULL = open(os.devnull, 'w')
#        qsub_return = subprocess.check_output(qsub_command,stderr=FNULL)
#        FNULL.close()
        qsub_return = subprocess.check_output(qsub_command,stderr=sys.stderr)
        print >>sys.stderr, "JOB NO: ",qsub_return.split('.')[0]
        #x = re.findall('Job <(\d+)> is submitted to queue <.*>',qsub_return)
        #no need if using -terse
#        x = qsub_return
#        if len(x)>0:
#            job_no = x[0]
        job_no=qsub_return.rstrip().split('.')[0]
    return (job_no, ' '.join(qsub_command))



#TO DO, WILL NOT WORK WITH GridEngine
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
    print >>killscript,"qdel ",job_id
    sublog.close()
    killscript.close()
    

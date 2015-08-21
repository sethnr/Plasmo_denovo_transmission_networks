#!/bin/bash

#SGE: run in same directory in which it is called
#$ -cwd

source /broad/software/scripts/useuse 

reuse -q Samtools
reuse -q BWA
reuse -q VCFtools
reuse -q Tabix
# use -q Python-3.4
reuse -q Python-2.7

reuse -q GCC-4.9
reuse -q Perl-5.10
reuse -q R-3.2
reuse -q Graphviz
reuse -q Discovar


COMMANDFILE=`readlink -f $1`
JOBINDEX=-1

#GET TAST ID
if [ -n "$LSB_JOBINDEX" ]; then
    JOBINDEX=${LSB_JOBINDEX}
elif [ -n "$SGE_TASK_ID" ]; then 
    JOBINDEX=${SGE_TASK_ID}
else
    echo "Need to set LSB_JOBINDEX or SGE_TASK_ID"
    exit 1
fi  

#GET JOB ID
if [ -n "$LSB_JOBID" ]; then
    JOBID=${LSB_JOBID}
elif [ -n "$JOB_ID" ]; then 
    JOBID=${JOB_ID}
else
    echo "Need to set LSB_JOBID or JOB_ID"
    exit 1
fi  

if [ -z $COMMANDFILE ]; then
    echo "Command file ${COMMANDFILE} not found"
    exit 1
fi  

COMMAND=`perl -ne "print if ${JOBINDEX}..${JOBINDEX}" $COMMANDFILE`

#LOG START OF JOB:
NOW=`date`
INSERT="db.dbjobs.insert({\"_id\":\"${JOBID}.${JOBINDEX}\", 
          \"JOB\":${JOBID},
          \"JOBINDEX\":${JOBINDEX},
          \"COMMAND\":\"${COMMAND}\",
          \"LAST\":new Date(),
          \"STATE\":\"STARTED\"});"
mongo $DBADD  --eval "${INSERT}"


#RUN COMMAND
echo "running command $JOBINDEX from file $COMMANDFILE"
echo $COMMAND;

${COMMAND}

JOBOUT=$?


#LOG COMPLETION BY DB:

if [ $JOBOUT == 0 ]; then
INSERT="db.dbjobs.update({\"_id\":\"${JOBID}.${JOBINDEX}\"},
          {\$set:{\"STATE\":\"DONE\"},
           \$currentDate:{\"LAST\":true}
          } );"
mongo $DBADD  --eval "${INSERT}"
else
INSERT="db.dbjobs.update({\"_id\":\"${JOBID}.${JOBINDEX}\"},
          {\$set:{\"STATE\":\"FAIL\"},
           \$currentDate:{\"LAST\":true}
          } );"
mongo $DBADD  --eval "${INSERT}"
fi

exit $JOBOUT

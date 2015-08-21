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

echo "running command $JOBINDEX from file $COMMANDFILE"
echo $COMMAND;

${COMMAND}

JOBOUT=$?

#LOG COMPLETION BY FILE:

#SWAP VARIABLE FOR COMMAND LOG (do edits elsewhere)
LOGFILE=${COMMANDFILE/txt/out.log}
if [[ ! -f $LOGFILE ]]; then
#cp $COMMANDFILE $LOGFILE
  touch $LOGFILE
fi

if [ $JOBOUT == 0 ]; then
   ( flock -x -w 120 200
   echo "#DONE:   ${COMMAND}" >>$LOGFILE
   ) 200>>$LOGFILE  || echo "DISCOVAR_HAS_NOT_LOGGED_PASS!" $?
else
   ( flock -x -w 120 200
   echo "#EXIT:   ${COMMAND}" >>$LOGFILE
   ) 200>>$LOGFILE  ||  echo "DISCOVAR_HAS_NOT_LOGGED_FAIL!" $?
fi

exit $JOBOUT

#!/bin/bash


COMMANDFILE=`readlink -f $1`
JOBINDEX=-1

if [ -n "$LSB_JOBINDEX" ]; then
    JOBINDEX=${LSB_JOBINDEX}
elif [ -n "$SGE_TASK_ID" ]; then 
    JOBINDEX=${SGE_TASK_ID}
else
    echo "Need to set LSB_JOBINDEX or SGE_TASK_ID"
    exit 1
fi  

if [ -z $COMMANDFILE ]; then
    echo "Command file ${COMMANDFILE} not found"
    exit 1
fi  


#COMMAND=`echo sed -n \'${LSB_JOB_INDEX}p\' < $COMMANDFILE`

COMMAND=`perl -ne "print if ${JOBINDEX}..${JOBINDEX}" $COMMANDFILE`

echo "running command $JOBINDEX from file $COMMANDFILE\n$COMMAND ";

${COMMAND}

DISCOOUT=$?

#SWAP VARIABLE FOR COMMAND LOG (do edits elsewhere)
# COMMANDFILE=${COMMANDFILE/txt/log}
LOGFILE=${COMMANDFILE/txt/out.log}


PRELEN=`wc -l $COMMANDFILE`
if [ $DISCOOUT == 0 ]; then

#    echo "DISCOVAR_HAS_WORKED!" $DISCOUT
#    echo cat "DONE\t ${COMMAND}" \>\>$LOGFILE
   ( flock -x -w 120 200
   echo "#DONE:   ${COMMAND}" >>$LOGFILE
   ) 200>>$LOGFILE  || echo "DISCOVAR_HAS_NOT_LOGGED_PASS!" $?

#   ( flock -x -w 120 200
#   perl -i -ne 'print "\#'${LSB_JOBINDEX}':DONE\t" if '${LSB_JOBINDEX}'..'${LSB_JOBINDEX}' ; print $_;' $COMMANDFILE
#   ) 200>>$COMMANDFILE  ||  echo "DISCOVAR_HAS_NOT_RECORDED_PASS!" $?
   
else
    
#    echo "DISCOVAR HAS FAILED!" $DISCOUT

   #( flock -x -w 120 200
   echo "#EXIT:   ${COMMAND}" >>$LOGFILE
   #) 200>>$LOGFILE  ||  echo "DISCOVAR_HAS_NOT_LOGGED_FAIL!" $?

#   ( flock -x -w 120 200
#   perl -i -ne 'print "\#'${LSB_JOBINDEX}':FAIL\t" if '${LSB_JOBINDEX}'..'${LSB_JOBINDEX}' ; print $_;' $COMMANDFILE
#   ) 200>>$COMMANDFILE  ||  echo "DISCOVAR_HAS_NOT_RECORDED_FAIL!" $?

fi

POSTLEN=`wc -l $COMMANDFILE`
if [[ $PRELEN != $POSTLEN ]]
then
    echo $COMMANDFILE " " $PRELEN " != " $POSTLEN
fi

exit $DISCOOUT

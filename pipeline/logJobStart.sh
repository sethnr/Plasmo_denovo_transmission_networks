#!/bin/bash

echo "logging job start"
JOBID=123456
JOBINDEX=2

DBADD='node1757.broadinstitute.org:27017/dbjobs'
NOW=`date`
COMMAND="did some stuff to some files"

INSERT="db.dbjobs.insert({\"_id\":\"${JOBID}.${JOBINDEX}\", 
          \"JOB\":${JOBID},
          \"JOBINDEX\":${JOBINDEX},
          \"COMMAND\":\"${COMMAND}\",
          \"LAST\":new Date(),
          \"STATE\":\"STARTED\"});"

echo $INSERT
echo mongo $DBADD  --eval "${INSERT}"
mongo $DBADD  --eval "${INSERT}"

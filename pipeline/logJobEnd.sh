#!/bin/bash

echo "logging job start"
JOBID=123456
JOBINDEX=2

DBADD='node1757.broadinstitute.org:27017/dbjobs'
NOW=`date`
INSERT="db.dbjobs.update({\"_id\":\"${JOBID}.${JOBINDEX}\"},
          {\$set:{\"STATE\":\"DONE\"},
           \$currentDate:{\"LAST\":true}
          } );"

echo $INSERT
echo mongo $DBADD  --eval "${INSERT}"
mongo $DBADD  --eval "${INSERT}"

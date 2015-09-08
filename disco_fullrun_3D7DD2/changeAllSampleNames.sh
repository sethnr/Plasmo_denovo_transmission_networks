#!/bin/bash

for MYS in `ls -d SM*1 SM*2`
do
  MYOLD=${MYS%_1}
  MYOLD=${MYOLD%_2}
  echo ../changeSampleNames.sh $MYS $MYOLD
done

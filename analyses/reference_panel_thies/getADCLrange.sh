#!/bin/bash

## defaults
TREE=""
CULTURED=""
NCULT=10
N=35
ITERATIONS=20
OUT="runLobSTR.out"

while getopts "t:C:c:n:o:q:" opt; do
  case $opt in
    t) TREE=${OPTARG} ;;
    C) CULTURED=$OPTARG ;;
    c) NCULT=$OPTARG ;;
    n) N=$OPTARG ;;
    i) ITERATIONS=$OPTARG;;
    o) OUT=$OPTARG ;;
    q) QUERY=$OPTARG;;
   \?) echo "Invalid option: -$OPTARG" >&2 ;;
    h) echo "  getADCLrange.sh -t <tree.newick> -C <CULTURED IDS>";
       echo "              -i <iterations> -c <no cultured to chose>";
       echo "              -n <no total leaves> -o <outfile>";

  esac
done


for i in $(seq 1 $ITERATIONS);
do
#    echo $i
    shuf -n $NCULT $CULTURED  | sort > ${CULTURED}_sub_${NCULT}_${i}

    rppr min_adcl_tree --leaves $N --algorithm pam \
	--always-include ${CULTURED}_sub_${NCULT}_${i} \
	--all-adcls-file ${CULTURED}_sub_${NCULT}_${i}.adcls \
        $TREE >  ${CULTURED}_sub_${NCULT}_${i}.removes
    grep -vf ${CULTURED}_sub_${NCULT}_${i}.removes  daniels_list.txt | sort > ${CULTURED}_sub_${NCULT}_${i}.keeps
    rm ${CULTURED}_sub_${NCULT}_${i}.removes
#        $TREE > /dev/null
    tail -n1 ${CULTURED}_sub_${NCULT}_${i}.adcls
done


shuf -n $NCULT $CULTURED  | sort > ${CULTURED}_sub_${NCULT}_0
for i in $(seq 1 $ITERATIONS);
do
#    echo $i

    rppr min_adcl_tree --leaves $N --algorithm pam \
	--always-include ${CULTURED}_sub_${NCULT}_0 \
	--all-adcls-file ${CULTURED}_sub_${NCULT}_0_${i}.adcls \
        $TREE >  ${CULTURED}_sub_${NCULT}_0_${i}.removes
    grep -vf ${CULTURED}_sub_${NCULT}_0_${i}.removes  daniels_list.txt | sort > ${CULTURED}_sub_${NCULT}_0_${i}.keeps
    rm ${CULTURED}_sub_${NCULT}_0_${i}.removes
#        $TREE > /dev/null
    tail -n1 ${CULTURED}_sub_${NCULT}_0_${i}.adcls
done


exit 0

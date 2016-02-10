#!/bin/bash

## defaults
TREE=""
CULTURED=""
NCULT=10
N=35
ITERATIONS=20
OUT="runLobSTR.out"

while getopts "t:C:c:n:o:q:i:" opt; do
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

CULTURED=`readlink -f $CULTURED`
TREE=`readlink -f $TREE`
TREENODES=`readlink -f daniels_list.txt`
mkdir tmp_${N}_sub_${NCULT}
cd tmp_${N}_sub_${NCULT}

echo $TREENODES
for i in $(seq 1 $ITERATIONS);
do
    echo $i
    shuf -n $NCULT $CULTURED  | sort > ${N}_sub_${NCULT}_${i}

    rppr min_adcl_tree --leaves $N --algorithm pam \
	--always-include ${N}_sub_${NCULT}_${i} \
	--all-adcls-file ${N}_sub_${NCULT}_${i}.adcls \
        $TREE >  ${N}_sub_${NCULT}_${i}.removes
    grep -vf ${N}_sub_${NCULT}_${i}.removes  $TREENODES | sort > ${N}_sub_${NCULT}_${i}.keeps
    rm ${N}_sub_${NCULT}_${i}.removes
#        $TREE > /dev/null
    tail -n1 ${N}_sub_${NCULT}_${i}.adcls
done
cat ${N}_sub_${NCULT}_*.adcls | perl -pe "s/\,/\t${NCULT}\t/gi" > ../${N}_sub_${NCULT}.adcls


#shuf -n $NCULT $CULTURED  | sort > ${CULTURED}_sub_${NCULT}_0
#for i in $(seq 1 $ITERATIONS);
#do
#    rppr min_adcl_tree --leaves $N --algorithm pam \
#	--always-include ${CULTURED}_sub_${NCULT}_0 \
#	--all-adcls-file ${CULTURED}_sub_${NCULT}_0_${i}.adcls \
#        $TREE >  ${CULTURED}_sub_${NCULT}_0_${i}.removes
#    grep -vf ${CULTURED}_sub_${NCULT}_0_${i}.removes  daniels_list.txt | sort > ${CULTURED}_sub_${NCULT}_0_${i}.keeps
#    rm ${CULTURED}_sub_${NCULT}_0_${i}.removes
#    tail -n1 ${CULTURED}_sub_${NCULT}_0_${i}.adcls
#done

exit 0

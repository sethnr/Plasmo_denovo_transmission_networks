#!/bin/bash

PED=$1

PEDTOTREE=$DISCO1/analyses/njtree_thies/getNJTreePedBOOT.R

Rscript $PEDTOTREE ${PED/.ped/}
if [[ $? != 0 ]]; then exit 1 ; fi


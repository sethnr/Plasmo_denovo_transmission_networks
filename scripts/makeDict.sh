#!/bin/bash

FASTA=$1
PICARD=/seq/software/picard/1.99/bin


java -jar ${PICARD}/CreateSequenceDictionary.jar R=${FASTA} O=${FASTA/.fasta/.dict}

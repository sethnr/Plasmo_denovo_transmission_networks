#!/bin/bash

export SV_DIR=~/software/svtoolkit/
SV_TMPDIR=./tmpdir

BAMLIST=$1 

INTERVALLIST=$WORK/refs/Pf3D7_v3.dict

mx="-Xmx4g"
BASEDIR=$PWD
WKDIR=${BAMLIST/.list/}

mkdir -p ${WKDIR}/logs || exit 1
mkdir -p ${WKDIR}/metadata || exit 1
cd ${WKDIR} || exit 1

classpath="${SV_DIR}/lib/SVToolkit.jar:${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar:${SV_DIR}/lib/gatk/Queue.jar"
java -Xmx4g -cp ${classpath} \
     org.broadinstitute.gatk.queue.QCommandLine \
     -S ${SV_DIR}/qscript/discovery/cnv/CNVDiscoveryPipeline.q \
     -S ${SV_DIR}/qscript/SVQScript.q \
     -gatk ${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar \
     -configFile ${SV_DIR}/conf/genstrip_parameters.txt \
     -R $WORK/refs/Pf3D7_v3.fasta \
     -I ${BASEDIR}/${BAMLIST} \
     -jobRunner Drmaa  -gatkJobRunner Drmaa \
     -md input_metadata_directory \
     -runDirectory ./ \
     -jobLogDir ./logs \
     -genderMapFile ${BASEDIR}/gender_map_file.txt \
     -ploidyMapFile ${BASEDIR}/ploidy_map_file.txt \
     -tilingWindowSize 1000 \
     -tilingWindowOverlap 500 \
     -maximumReferenceGapLength 1000 \
     -boundaryPrecision 100 \
     -minimumRefinedLength 500 \
     -run

#     -intervalList reference_chromosomes.list \

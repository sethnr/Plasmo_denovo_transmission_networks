#!/bin/bash

for VCF in $@
do
    vcf-sort -c $VCF > ${VCF}.tmp
    if [[ $? != 0 ]]; then exit $?; fi
    mv ${VCF}.tmp ${VCF}
    bgzip ${VCF}
    if [[ $? != 0 ]]; then exit $?; fi
    tabix -pvcf ${VCF}.gz
done

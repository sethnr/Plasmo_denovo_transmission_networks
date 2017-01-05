#!/bin/bash 

TREEFOLDER=/seq/plasmodium/sredmond/pfdisco/analyses/njtree_thies/

for SET in ALL SNP INDEL CALLBOTH CALLDISCOONLY
do
    for SUMM in mcct
    do
	python sumtrees.py -F nexus -s $SUMM \
	    --output-tree-filepath=${TREEFOLDER}/boot_${SET}/sum_${SET}.${SUMM}.nexus \
	    ${TREEFOLDER}/boot_${SET}/Thies*nexus
    done
    
#	-t ${TREEFOLDER}/Thies_all_manual.PASS.Cls.miss0.5.LMRG.HAP.vcf.NJ.nexus
#	-t ${TREEFOLDER}/Thies_manual.nj.NJ.nexus \
    python sumtrees.py -F nexus \
	--output-tree-filepath=${TREEFOLDER}/boot_${SET}/sum_${SET}.target.nexus \
	-t ${TREEFOLDER}/Thies_all_manual.PASS.Cls.miss0.5.LMRG.HAP.vcf.NJ.nexus \
	${TREEFOLDER}/boot_${SET}/Thies*nexus
   
done



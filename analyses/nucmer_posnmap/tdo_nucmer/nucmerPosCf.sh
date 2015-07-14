#!/bin/bash

name=$1
ref=$2
qry=$3

nucmer -p $name $ref $qry
delta-filter -m -i 97 -l 1000 $name.delta >$name.filter.delta
show-coords -THqcl -o $name.filter.delta> $name.filter.coords
show-snps -CHTr $name.filter.delta > $name.Csnp
show-snps -HTr $name.filter.delta > $name.snp 

y=$name;
perl ~tdo/Bin/reichenowi.allels.V2.pl $y.filter.coords $y.Csnp $y.snp $ref3D7 $y  > Nuc.alleles.$name.txt


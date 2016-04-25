#!/bin/bash

for V1 in disco haplo100 haplo250 ;
do
  for V2 in disco haplo100 haplo250;
    do
      vcftools --not-chr M76611 --vcf ${V1}.M50.vcf --diff ${V2}.M50.vcf --diff-site --out ${V1}${V2}.M50
      vcftools --not-chr M76611 --vcf ${V1}.vcf --diff ${V2}.vcf --diff-site --out ${V1}${V2}
  done
done

#!/bin/bash

pheno=$1
mkdir $pheno

for i in {1..22}; do 
    ./plink --bfile /mnt/project/Bulk/'Genotype Results'/'Genotype calls'/ukb22418_c${i}_b0_v2 \
    --clump-p1 1e-7 --clump-p2 1e-2 --clump-r2 0.2 --clump-kb 250 \
    --clump common.assoc.tsv \
    --out ./$pheno/c${i} 
done

awk -F" " '{print $1,$4,$3}' ./$pheno/*.clumped > $pheno.common.loci.txt

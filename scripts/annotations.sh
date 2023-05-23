#!/bin/bash

# Loop over chromosomes 1 to 22
for chr in {11..22}; do
    echo "Processing chromosome ${chr}..."
    
    # Run LD score regression for the current chromosome
    python ldsc/ldsc.py \
        --l2 \
        --bfile data/1000g_EUR_Phase3_plink/1000G.EUR.QC.${chr} \
        --ld-wind-cm 1 \
        --annot data/bscores/snpwise_bim/snpwise.CADD.chr${chr}.bmap.txt \
        --out data/bscores/ld_score/snpwise.CADD.chr${chr}
done


python ldsc/ldsc.py \
--h2 data/GWAS_summaries/munged/munged.4125.info0.allchroms.pan.sumstats.gz \
--ref-ld data/bscores/ld_score/snpwise.CADD.chr3 \
--w-ld-chr data/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
--frqfile data/1000G_Phase3 \
--print-coefficients \
--print-delete-vals \
--print-cov \
--return-silly-things \
--out data/results/rg/tSDS.4125.info0.allchroms.pan \
--M 7009236 \
--overlap-annot 

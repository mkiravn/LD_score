##### A snakemake pipeline to run LD score regression

#wildcard_constraints:
#   pheno_code="^\d+(?:_[a-z]*)*$"



rule all:
    input:
        expand("results/{pheno_code}.{spop}.info{info}.chr{chrom}.log",info=0,spop="EUR",chrom=1,pheno_code=[100890,3731,"50_raw"])


#This rule unzips raw bgz files
rule unzip_bgz:
   input:
       "data/GWAS_summaries/raw/{file}.tsv.bgz"
   output:
       "data/GWAS_summaries/unzipped/{file}.tsv"
   shell:
       "gunzip -c {input} > {output}"

rule download_sumstats_files:
    """
    Download phenotype files
    """
    output:
        "data/GWAS_summaries/raw/{pheno_code}.gwas.imputed_v3.both_sexes.tsv.bgz"
    shell:
       "wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/{wildcards.pheno_code}.gwas.imputed_v3.both_sexes.tsv.bgz -O data/GWAS_summaries/raw/{wildcards.pheno_code}.gwas.imputed_v3.both_sexes.tsv.bgz"

rule download_ldscore_files:
    """
    Download ld score files
    """
    output:
        "data/UKBB.ALL.ldscore.tar.gz" # sub-optimal, should be all pops we want
    shell:
        "wget https://pan-ukb-us-east-1.s3.amazonaws.com/ld_release/UKBB.ALL.ldscore.tar.gz > {output}"

rule unzip_ldscore_files:
	input:
		"data/UKBB.ALL.ldscore.tar.gz"
	output: 
		"data/UKBB.ALL.ldscore/UKBB.EUR.rsid.l2.ldscore.gz"
	shell: 
		"tar –xvzf {input}"

rule download_variant_file:
    """
    Download variant file
    """
    output:
        "data/GWAS_summaries/raw/variants.tsv.bgz"
    shell:
       "wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/annotations/variants.tsv.bgz -O data/GWAS_summaries/raw/variants.tsv.bgz"


rule filter_variants:
    """
    Filter variants by INFO score and chromosome
    """
    input:
        "data/GWAS_summaries/unzipped/variants.tsv"
    output:
        "data/GWAS_summaries/analysis/variants.info{info}.chr{chrom}.tsv"
    shell:
        "Rscript scripts/filter_variants.R {input} {output} {wildcards.info} {wildcards.chrom}"


rule combine_sumstats:
    """
    Combine info from variant file with info from summary statistics file
    """
    input:
        variants="data/GWAS_summaries/analysis/variants.info{info}.chr{chrom}.tsv",
        sumstats= "data/GWAS_summaries/unzipped/{pheno_code}.gwas.imputed_v3.both_sexes.tsv"
    output:
        combined="data/GWAS_summaries/unzipped/{pheno_code}.info{info}.chr{chrom}.sumstats.tsv"
    shell:
        "Rscript scripts/combine_sumstats.R {input.variants} {input.sumstats} {output.combined}"

rule munge_sumstats:
    """
    Processing the summary statistics file into a format which ldsc likes
    """
    input:
        "data/GWAS_summaries/unzipped/{pheno_code}.info{info}.chr{chrom}.sumstats.tsv"
    output:
        "data/GWAS_summaries/sumstats/munged.{pheno_code}.info{info}.chr{chrom}.sumstats.tsv"
    params:
        output_path = "data/GWAS_summaries/sumstats/munged.{pheno_code}.info{info}.chr{chrom}"
    conda:
        "ldsc/environment.yml"
    shell: # this is a *very* hacky way to get it to activate conda...
        """
        source activate ldsc && \
        ./ldsc/munge_sumstats.py \
        --out {params.output_path} \
        --sumstats {input} \
        --N 361194.0 --a1 alt --a2 ref
        """

# ./ldsc/munge_sumstats.py --out data/GWAS_summaries/sumstats/munged.100890.info0.chr1 --sumstats data/GWAS_summaries/sumstats/100890.info0.chr1.sumstats.tsv --N 361194.0 --a1 alt --a2 ref
# worked on the command line???

rule cut_ld_scores:
    """Format LD files into by cutting only the baseL2 column"""
    input:
        ld_score = "data/UKBB.ALL.ldscore/UKBB.{spop}.8LDMS.rsid.l2.ldscore.gz"
    output:
        "data/ld_score/cut.UKBB.{spop}.8LDMS.rsid.l2.ldscore"
    params:
         output_dir = "data/ld_score/"
    shell:    
        """
        mkdir -p {output_dir} && \
        gunzip -c {input} | \
        cut -f 1-4 > {output}
        """


rule estimate_heritability:
    """
    Using LD score regression to estimate heritability
    """
    input:
        sumstats = "data/GWAS_summaries/sumstats/munged.{pheno_code}.info{info}.chr{chrom}.sumstats.tsv",
        ld = "data/UKBB.ALL.ldscore/UKBB.{spop}.rsid.l2.ldscore.gz"
    output:
        "results/{pheno_code}.{spop}.info{info}.chr{chrom}.log"
    params:
        ld_path = "data/UKBB.ALL.ldscore/UKBB.{spop}",
        output_path = "results/{pheno_code}.info{info}.chr{chrom}",
        M = 7009236 # modify this
    conda:
        "ldsc/environment.yml"
    shell:
        """
        ./ldsc/ldsc.py \
        --h2 {input.sumstats} \
        --ref-ld {params.ld_path} \
        --w-ld {params.ld_path} \
        --M {params.M} \
        --out {params.output_path}
        """

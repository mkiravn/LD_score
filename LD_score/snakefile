##### A snakemake pipeline to run LD score regression


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
        "data/GWAS_summaries/sumstats/{pheno_code}.info{info}.chr{chrom}.sumstats.tsv"
    shell:
        "Rscript scripts/combine_sumstats.R {input.variants} {input.sumstats} {output}"

rule munge_sumstats:
    """
    Processing the summary statistics file into a format which ldsc likes
    """
    input:
         "data/GWAS_summaries/sumstats/{pheno_code}.info{info}.chr{chrom}.sumstats.tsv"
    output:
        "data/GWAS_summaries/sumstats_munged/{pheno_code}.info{info}.chr{chrom}.sumstats.tsv"
    params:
        output_path= "data/GWAS_summaries/sumstats_munged/{pheno_code}.info{info}.chr{chrom}"
    conda:
        "ldsc/environment.yml"
    shell:
        """./ldsc/munge_sumstats.py \
                --out {params.output_path} \
                --N 361194.0 \
                --sumstats {input} """
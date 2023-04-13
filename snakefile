##### A snakemake pipeline to run LD score regression

#wildcard_constraints:
#   pheno_code="^\d+(?:_[a-z]*)*$"

import pandas as pd

# Read the TSV file into a pandas DataFrame
df = pd.read_csv("data/simons_phenotypes.tsv", sep="\t")

# Extract the values of the "code" column into a list
code_list = df["code"].tolist()[:1]
ns = df["N"].tolist()[:1]
# note to self - deal with N (in phenotype file, currently ignored)

rule all:
    """will return the input files needed for ldsc"""
    input:
        expand("data/GWAS_summaries/processed/{pheno_code}.info{info}.allchroms.sumstats.tsv",info=0,pheno_code=code_list),
        expand("data/sds/{pheno_code}.info{info}.allchroms.tSDS.tsv",info=0,pheno_code=code_list)

# rule all:
#     """goes up until the end of the wgets"""
#     input:
#         "data/UKBB.ALL.ldscore.tar.gz",
#         "data/sds/SDS_UK10K_n3195_release_Sep_19_2016.tab.gz",
#         expand("data/GWAS_summaries/raw/{pheno_code}.gwas.imputed_v3.both_sexes.tsv.bgz",pheno_code=code_list),
        


rule download_sds:
    output:
        "data/sds/SDS_UK10K_n3195_release_Sep_19_2016.tab.gz"
    shell:
        "wget https://datadryad.org/stash/downloads/file_stream/10751 -O data/sds/SDS_UK10K_n3195_release_Sep_19_2016.tab.gz"

rule unzip_sds:
    input:
        "data/sds/SDS_UK10K_n3195_release_Sep_19_2016.tab.gz"
    output:
        "data/sds/SDS_UK10K_n3195_release_Sep_19_2016.tab"
    shell:
        "gunzip -c {input} > {output}"

rule download_bscores:
    output:
        "data/bscores/CADD_bestfit.tar.gz"
    shell:
        "wget https://github.com/sellalab/HumanLinkedSelectionMaps/raw/master/Bmaps/CADD_bestfit.tar.gz -O data/bscores/CADD_bestfit.tar.gz"


# these rules only get the phenotype table - ironically though these are needed in my current rule all. 
rule download_phenotype_table:
    """downloads table of phenotypes from Simons et al. 2023, 
    https://www.biorxiv.org/content/10.1101/2022.10.04.509926v1.supplementary-material"""
    output:
        "data/simons_phenotypes.csv"
    shell:
        "wget https://www.biorxiv.org/content/biorxiv/early/2022/10/07/2022.10.04.509926/DC2/embed/media-2.csv?download=true -O data/simons_phenotypes.csv"

rule process_phenotype_table:
    input:
       "data/simons_phenotypes.csv" 
    output:
       "data/simons_phenotypes.tsv" 
    shell:
        "Rscript scripts/get_phenotypes.R {input} {output}"

# rule unzip_bscores:
#     input:
#         "data/bscores/CADD_bestfit.tar.gz"
# 	output: 
#         "data/bscores/CADD_bestfit/ch1.bmap.txt"
#     params:
#         path="data/bscores/CADD_bestfit"
# 	shell: 
#         "tar -xzf {input} -O {params.path}"


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
       """
       mkdir -p data/GWAS_summaries/raw/ && \
       wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/{wildcards.pheno_code}.gwas.imputed_v3.both_sexes.tsv.bgz \
       -O data/GWAS_summaries/raw/{wildcards.pheno_code}.gwas.imputed_v3.both_sexes.tsv.bgz
       """

rule download_infant_mortality_sumstats:
    """rule to download summary statistics for infant mortality"""
    output:
        "data/GWAS_summaries/raw/Sumstats_byimr_Wu_2021Jun.txt.gz"
    shell:
        """
        wget ftp://ftp.biostat.wisc.edu/pub/lu_group/Projects/IMR/Sumstats/Sumstats_byimr_Wu_2021Jun.txt.gz \
        -O data/GWAS_summaries/raw/Sumstats_byimr_Wu_2021Jun.txt.gz &&
        gunzip -c data/GWAS_summaries/raw/Sumstats_byimr_Wu_2021Jun.txt.gz 
        """


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
       """
       mkdir -p data/GWAS_summaries/raw/ && \
       wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/annotations/variants.tsv.bgz -O \
       data/GWAS_summaries/raw/variants.tsv.bgz
       """


rule filter_variants:
    """
    Filter variants by INFO score and chromosome
    """
    input:
        "data/GWAS_summaries/unzipped/variants.tsv"
    output:
        "data/GWAS_summaries/analysis/variants.info{info}.chr{chrom}.tsv"
    shell:
        r"""
        awk -F '\t' 'BEGIN {{OFS="\\t"}} NR==1 {{print}} ($2 == {wildcards.chrom}) && ($7 > {wildcards.info}) {{print}}' \
        data/GWAS_summaries/unzipped/variants.tsv > "data/GWAS_summaries/analysis/variants.info{wildcards.info}.chr{wildcards.chrom}.tsv"
        """
#        "Rscript scripts/filter_variants.R {input} {output} {wildcards.info} {wildcards.chrom}"


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

rule polarise_sds:
    """polarises and munges sds"""
    input:
        sds="data/sds/SDS_UK10K_n3195_release_Sep_19_2016.tab",
        sumstats="data/GWAS_summaries/unzipped/{pheno_code}.info{info}.chr{chrom}.sumstats.tsv",
    output:
        polarised_sds="data/sds/{pheno_code}.info{info}.chr{chrom}.tSDS.tsv"
    shell:
        """
        Rscript scripts/polarise_SDS.R {input.sds} {input.sumstats} {output} 
        """ 

rule pool_sumstats:
    """pools summary statistics into one file"""
    input:
        expand("data/GWAS_summaries/unzipped/{pheno_code}.info{info}.chr{chrom}.sumstats.tsv", pheno_code="{pheno_code}", info="{info}", chrom=range(1, 23))
    output:
        "data/GWAS_summaries/processed/{pheno_code}.info{info}.allchroms.sumstats.tsv"
    shell:
        """
        find data/GWAS_summaries/unzipped/ -name "{wildcards.pheno_code}.info{wildcards.info}.chr*.sumstats.tsv" -type f -print0 | \
        xargs -0 awk -F'\\t' 'FNR==1 && NR>1 {{next}} 1' | \
        sort -n -k4,4 -k5,5 > {output}
        """

rule pool_SDS:
    """pools summary SDS into one file"""
    input:
        expand("data/sds/{pheno_code}.info{info}.chr{chrom}.tSDS.tsv", pheno_code="{pheno_code}", info="{info}", chrom=range(1, 23))
    output:
        "data/sds/{pheno_code}.info{info}.allchroms.tSDS.tsv"
    shell:
        """
        find data/sds -name "{wildcards.pheno_code}.info{wildcards.info}.chr*.tSDS.tsv" -type f -print0 | \
        xargs -0 awk 'FNR==1 && NR>1 {{next}} 1'| \
        sort -n -k1,1 -k2,2 > {output}
        """

# # works up to here...
# rule munge_sumstats:
#     """
#     Processing the summary statistics file into a format which ldsc likes
#     """
#     input:
#         "data/GWAS_summaries/unzipped/{pheno_code}.info{info}.chr{chrom}.sumstats.tsv"
#     output:
#         "data/GWAS_summaries/sumstats/munged.{pheno_code}.info{info}.chr{chrom}.sumstats.tsv"
#     params:
#         output_path = "data/GWAS_summaries/sumstats/munged.{pheno_code}.info{info}.chr{chrom}"
#     conda:
#         "ldsc/environment copy.yml"
#     shell: # this is a *very* hacky way to get it to activate conda...
#         """
#         ./ldsc/munge_sumstats.py \
#         --out {params.output_path} \
#         --sumstats {input} \
#         --N 361194.0 --a1 alt --a2 ref &&  
#         """


# rule estimate_heritability:
#     """
#     Using LD score regression to estimate heritability - currently does not work
#     """
#     input:
#         sumstats = "data/GWAS_summaries/sumstats/munged.{pheno_code}.info{info}.chr{chrom}.sumstats.tsv",
#         ld = "data/UKBB.ALL.ldscore/UKBB.{spop}.rsid.l2.ldscore.gz"
#     output:
#         "results/{pheno_code}.{spop}.info{info}.chr{chrom}.log"
#     params:
#         ld_path = "data/UKBB.ALL.ldscore/UKBB.{spop}",
#         output_path = "results/{pheno_code}.info{info}.chr{chrom}",
#         M = 7009236 # modify this
#     conda:
#         "ldsc/environment.yml"
#     shell:
#         """
#         ./ldsc/ldsc.py \
#         --h2 {input.sumstats} \
#         --ref-ld {params.ld_path} \
#         --w-ld {params.ld_path} \
#         --M {params.M} \
#         --out {params.output_path}
#         """

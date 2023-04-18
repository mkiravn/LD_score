#!/bin/bash
# read in codes
codes=$(cut -f 2 data/simons_phenotypes.tsv | head -n 2 | awk 'NR > 1')
# define input files
input_files=(
    "50_raw.info0.allchroms"
)

modified_codes=""
for code in $codes
do
    modified_code="$code.info0.allchroms"
    # Append modified code to variable
    modified_codes="$modified_codes $modified_code"
done

# Print modified codes to console
echo "Running analysis for phenotypes: $modified_codes"

# define directories and files
input_dir="data/GWAS_summaries/processed"
munged_dir="data/GWAS_summaries/munged"
output_dir="data/results/rg"
sds_dir="data/sds"

# loop over input files and execute scripts
for input_file in ${modified_codes}; do
    # generate munged file path
    munged_file="${munged_dir}/munged.${input_file}"
    # generate output file path
    output_file="${output_dir}/${input_file}"
    sds_file="${sds_dir}/${input_file}.tSDS.tsv"

    # run the scripts
    mkdir -p ${output_dir}
    mkdir -p ${munged_dir}
    bash scripts/munge_stuff.sh "${input_dir}/${input_file}.sumstats.tsv" "${munged_file}" && \
    bash scripts/estimate_rg.sh "${munged_file}.sumstats.gz" "${sds_file}" data/UKBB.ALL.ldscore/UKBB.EUR.rsid "${output_file}" 23960350
done

#!/bin/bash
# read in codes
codes=$(cut -f 2 data/simons_phenotypes.tsv | awk 'NR > 1')
# define input files
input_files=(
    "50_raw.info0.chr1"
    "3731.info0.chr1"
    "100890.info0.chr1"
)

modified_codes=""
for code in $codes
do
    modified_code="$code.info0.chr1"
    # Append modified code to variable
    modified_codes="$modified_codes $modified_code"
done

# Print modified codes to console
echo "Running analysis for phenotypes: $modified_codes"

# define directories and files
input_dir="data/GWAS_summaries/unzipped"
munged_dir="data/GWAS_summaries/sumstats"
output_dir="data/results/rg"

# loop over input files and execute scripts
for input_file in ${modified_codes}; do
    # generate munged file path
    munged_file="${munged_dir}/munged.${input_file}"
    # generate output file path
    output_file="${output_dir}/${input_file}"
    
    # run the scripts
    mkdir -p data/results
    bash scripts/munge_stuff.sh "${input_dir}/${input_file}.sumstats.tsv" "${munged_file}" && \
    bash scripts/estimate_rg.sh "${munged_file}.sumstats.gz" "${input_file}.tSDS.tsv" EUR "${output_file}"
done

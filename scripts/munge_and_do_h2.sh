#!/bin/bash

# define input files
input_files=(
    "50_raw.info0.chr1"
    "3731.info0.chr1"
    "100890.info0.chr1"
)

# define directories and files
input_dir="data/GWAS_summaries/unzipped/"
munged_dir="data/GWAS_summaries/sumstats"
output_dir="data/results"

# loop over input files and execute scripts
for input_file in "${input_files[@]}"; do
    # generate munged file path
    munged_file="${munged_dir}/munged.${input_file}"
    # generate output file path
    output_file="${output_dir}/${input_file}"
    
    # run the scripts
    mkdir -p data/results
    bash scripts/munge_stuff.sh "${input_dir}/${input_file}.sumstats.tsv" "${munged_file}" && \
    bash scripts/estimate_h2.sh "${munged_file}.sumstats.gz" EUR "${output_file}"
done

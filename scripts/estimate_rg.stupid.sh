#!/bin/bash

# Get the input file from the command line
sumstats_file=$1
tsds_file=$2
ld_path=$3
output_path=$4
M=$5




# Run the ldsc command with input file and other parameters
./ldsc/ldsc.py \
--rg "${sumstats_file}","${tsds_file}" \
--ref-ld "${ld_path}" \
--w-ld "${ld_path}" \
--M 10000 \
--out "${output_path}" \
--return-silly-things \
--print-cov


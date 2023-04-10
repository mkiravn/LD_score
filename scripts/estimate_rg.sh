#!/bin/bash

# Get the input file from the command line
sumstats_file=$1
tsds_file=$3
spop=$4
output_path=$4

# Set other parameters
M=361194.0


# Run the ldsc command with input file and other parameters
./ldsc/ldsc.py \
--rg "${sumstats_file}","${tsds_file}" \
--ref-ld "data/UKBB.ALL.ldscore/UKBB.${spop}.rsid" \
--w-ld "data/UKBB.ALL.ldscore/UKBB.${spop}.rsid" \
--M "${M}" \
--out "${output_path}"

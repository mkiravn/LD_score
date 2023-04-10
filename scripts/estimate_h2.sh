#!/bin/bash

# Get the input file from the command line
input_file=$1
spop=$2
output_path=$3

# Set other parameters
M=361194.0


# Run the ldsc command with input file and other parameters
./ldsc/ldsc.py \
--h2 "${input_file}" \
--ref-ld "data/UKBB.ALL.ldscore/UKBB.${spop}.rsid" \
--w-ld "data/UKBB.ALL.ldscore/UKBB.${spop}.rsid" \
--M "${M}" \
--out "${output_path}"

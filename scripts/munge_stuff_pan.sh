#!/bin/bash

# Get the input and output files from the command line
input_file=$1
munged_file=$2

# Set other parameters
N=361194.0
a1="A1"
a2="A2"

# Run the munge_sumstats.py command with input and output files
./ldsc/munge_sumstats.py \
--out "${munged_file}" \
--sumstats "${input_file}" \
--N-col "n_cases_EUR" \
--a1 "${a1}" \
--a2 "${a2}" \
--p "neglog10_pval_EUR" \
--signed-sumstats beta_EUR,0

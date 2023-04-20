library(tidyverse)
library(data.table)


# Read in command line arguments
args <- commandArgs(trailingOnly=TRUE)

# read in the input files
variants <- args[1]
sumstats <- args[2]
phenotypes <- args[3]
code <- args[4]

print(paste("Reading in files:",variants,sumstats))

n <- phenotypes %>% filter(phenocode==code)

variants <- fread(variants,header=TRUE) %>% rename("chr"="chrom") %>% mutate(chr=as.integer(chr))
sumstats <- fread(sumstats,header=TRUE) %>% mutate(chr=as.integer(chr))

# output path name
output_path <- args[5]


# join and write out
print(paste("Processing. Will output to:",output_path))
left_join(variants,sumstats,by=c("chr","pos","ref","alt")) %>% # different with the other dataset
    mutate(A1=alt,A2=ref,n_cases_EUR=n$n_cases_EUR) %>%
    select(chr,pos,ref,alt,A1,A2,rsid,varid,beta_EUR,neglog10_pval_EUR,n_cases_EUR) %>%
    fwrite(sep="\t", file=output_path)



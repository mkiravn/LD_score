library(tidyverse)
library(data.table)


# Read in command line arguments
args <- commandArgs(trailingOnly=TRUE)

# read in the input files
variants <- args[1]
sumstats <- args[2]

print(paste("Reading in files:",variants,sumstats))

variants <- fread(variants,header=TRUE)
sumstats <- fread(sumstats,header=TRUE)

# output path name
output_path <- args[3]


# join and write out
print(paste("Processing. Will output to:",output_path))
left_join(variants,sumstats,by=c("chr","pos")) %>%
    mutate(a1=alt,a2=ref) %>%
    select(ref,alt,rsid,varid,chr,pos,beta,pval,n_complete_samples) %>%
    fwrite(sep="\t", file=output_path)



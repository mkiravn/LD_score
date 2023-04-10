library(tidyverse)
library(data.table)
### does some very low key reformatting

# Read in command line arguments
args <- commandArgs(trailingOnly=TRUE)

# read in the input files
phenotypes <- args[1]
output_path <- args[2]

phenos <- read.csv(phenotypes,header=TRUE)

phenos <- phenos %>% select(Trait,code=`UKBB.code`,N=`Median.study.size`,h2=`Estimated.heritability..h.2.`)

phenos %>% write_delim(file=output_path,delim="\t")

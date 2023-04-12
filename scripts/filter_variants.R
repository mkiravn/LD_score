library(tidyverse)
library(data.table)


# Read in command line arguments
args <- commandArgs(trailingOnly=TRUE)

# Assign the first argument to the input path variable
input_path <- args[1]

# Assign the second argument to the output path variable
output_path <- args[2]

# minimum info score
info_score <- args[3]

# chromosome 
chrom <- args[4]


# Read in the input file
input_data <- fread(input_path, header=TRUE)
print(paste("Read in input:",input_path))

# Filter by info score and chromosome if provided
if (!chrom=="all") {
    input_data <- input_data %>% filter(info > info_score, chr == chrom)
    print(paste("Filtered by chromosome. Writing out file at:", output_path))
} else {
    input_data <- input_data %>% filter(info>info_score)
    print(paste("Writing out file at:", output_path))
}


print(paste("Filtered. Writing out file at:",output_path))
input_data %>% fwrite(sep="\t", file=output_path)

library(tidyverse)
# Directory containing the files
directory <- 'data/results/rg'

## read in phenotypes - will allow us to read in codes 
# phenotypes <- read.table("data/simons_phenotypes.tsv", header = TRUE, sep = "\t")

# Name pattern of the files
filename_pattern <- 'tSDS.*.log'

output_file <-'data/results/rg/combined_rg.tsv'

# List all files in the directory
files <- list.files(directory, pattern = filename_pattern, full.names = TRUE)
print(paste("Files:",files))

# Initialize an empty data frame to store the combined results
combined_results <- data.frame()

# Loop through the files and read the results
for (file in files) {

    lines <- readLines(file)
  
    # Extract the lines containing the data rows
    data_lines <- lines[65]
    results <- data.frame(do.call(rbind, strsplit(data_lines, '\\s+')))
    results <- results[c(2:length(results))]
    # Append the results to the combined results data frame
    combined_results <- rbind(combined_results, results)
  
    # Extract the lines containing the data rows
    head_lines <- lines[64]
    header <- strsplit(head_lines, '\\s+')[[1]]
}
colnames(combined_results) <- header[c(2:length(header))]
# Print the combined results
write_tsv(combined_results,file=output_file)

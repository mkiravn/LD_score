library(tidyverse)
# Directory containing the files
directory <- 'data/results/rg'

## read in phenotypes - will allow us to read in codes 
# phenotypes <- read.table("data/simons_phenotypes.tsv", header = TRUE, sep = "\t")

# Name pattern of the files
filename_pattern <- 'tSDS.*.info0.allchroms.pan.log'

output_file <-'data/results/rg/combined_rg_pan.tsv'

# List all files in the directory
files <- list.files(directory, pattern = filename_pattern, full.names = TRUE)
print(paste("Files:",files))

# Initialize an empty data frame to store the combined results
combined_results <- data.frame()

# Loop through the files and read the results
for (file in files) {
    
        # Read the file into a character vector
	file_lines <- readLines(file)

	# Find the line number of "Summary of Genetic Correlation Results"
	summary_line <- grep("Summary of Genetic Correlation Results", file_lines)

	# Extract the two lines underneath "Summary of Genetic Correlation Results"
	results_lines <- file_lines[(summary_line+2)]
        
	# Split the lines by tab to create a data frame
	results_df <- t(data.frame(strsplit(results_lines, "\\s+"), stringsAsFactors = FALSE))
	
	colnames(results_df) <- c("p1", "p2", "rg", "se", "z", "p", "h2_obs", "h2_obs_se", "h2_int", "h2_int_se", "gcov_int", "gcov_int_se")

	# Extract the first line as column names
	#colnames(results_df) <- strsplit(file_lines[(summary_line+1)], "\t")[[1]]
	#results_df <- data.frame(read.table(results_lines,header=T))
	print(results_df)
  
    # Append the results to the combined results data frame
    combined_results <- rbind(combined_results, results_df)
  
}
# Print the combined results
write_tsv(combined_results,output_file)

options(warn = 0)
library(tidyverse)
# Read in codes
phenotypes <- read.table("data/simons_phenotypes_codes.tsv", header = TRUE, sep = "\t") %>% head(20)

codes <- phenotypes$ncode

# Define directories and files
input_dir <- "data/GWAS_summaries/processed"
munged_dir <- "data/GWAS_summaries/munged"
tSDS_dir <- "data/sds"
output_dir <- "data/results/rg"
ld_path <- "data/UKBB.ALL.ldscore/UKBB.EUR.rsid"
if (!file.exists(paste0(ld_path,".l2.ldscore.gz"))) {
    print(paste("LD file not found:", paste0(ld_path,".l2.ldscore.gz")))
    stop 
}


dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Loop over input files and execute scripts
for (i in c(1:dim(phenotypes)[1])) {
    code <- phenotypes$ncode[i]
    #print(paste("Running for code:",code))

    modified_code <- paste0(code, ".info0.allchroms.pan") # Modify codes
    # Generate input file paths
    sumstats_file <- file.path(input_dir, paste0(modified_code, ".sumstats.tsv"))
    tSDS_file <- file.path(tSDS_dir,paste0(modified_code, ".tSDS.tsv"))
    # Generate output file paths
    output_file <- file.path(output_dir, paste0("tSDS.",modified_code))
    munged_file_prefix <- file.path(munged_dir, paste0("munged.", modified_code))
    
    # Assert all the input files exist:

    # Check if sumstats file exists
    if (!file.exists(sumstats_file)) {
        print(paste("Sumstats file not found:", sumstats_file))
        next  # Skip to the next iteration of the loop
    }

    # Check if tSDS file exists
    if (!file.exists(tSDS_file)) {
        print(paste("tSDS file not found:", tSDS_file))
        next  # Skip to the next iteration of the loop
    }


    # Run the scripts
    system(paste0("bash scripts/munge_stuff_pan.sh ", sumstats_file, " ", munged_file_prefix))
    system(paste("bash scripts/estimate_rg.sh", paste0(munged_file_prefix, ".sumstats.gz"), tSDS_file, ld_path, output_file, 7009236))
}


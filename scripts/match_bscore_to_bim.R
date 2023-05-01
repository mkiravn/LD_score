###### Converts B-score format to SNP-wise scores by inner join on bim file
###### 
###### 

library(tidyverse)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

dir.create(file.path("data/bscores/snpwise_bim/"))
input_file <- args[1]
bim_file <- args[2]
output_file <- args[3]


bmap <- fread(input_file)
bim <- fread(bim_file,  header=F)
colnames(bim) <- c("CHR", "SNP","X","BP","ref","alt")

colnames(bmap) <- c("B","length")

bmap <- bmap %>%
          mutate(bin_end=cumsum(length),bin_start=bin_end-length,B=B/1000,CHR=chr)

# Define function to match SNP BP to bin in bmap
match_snp_to_bin <- function(snp_bp, bmap) {
  i <- which.max(bmap$bin_end > snp_bp)
  if (i == 0 || is.na(bmap$bin_start[i]) || is.na(bmap$bin_end[i])) {
    return(NA)
  } else {
    return(bmap$bin_end[i])
  }
}

bim <- bim %>% filter(CHR==chr) %>%
  rowwise() %>%
  mutate(bin_end=match_snp_to_bin(BP,bmap=bmap))

# Merge ldscore and bmap by bin
merged <- bim %>% left_join(bmap,by=c("bin_end","CHR"))

merged %>% select(CHR,SNP,BP,B) %>% write_tsv(file=output_file)

###### Converts B-score format to SNP-wise scores
###### 
###### 

library(tidyverse)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

dir.create(file.path("data/bscores/snpwise/"))
input_file <- args[1]
ld_file <- args[2]
output_path <- args[3]
chr <- args[4]

bmap <- fread(input_file)
ldscore <- fread(ld_file)

colnames(bmap) <- c("B","length")

bmap <- bmap %>%
          mutate(bin_end=cumsum(length),bin_start=bin_end-length,B=B/1000,CHR=as.numeric(chr))

# Define function to match SNP BP to bin in bmap
match_snp_to_bin <- function(snp_bp, bmap) {
  i <- which.max(bmap$bin_end > snp_bp)
  if (i == 0 || is.na(bmap$bin_start[i]) || is.na(bmap$bin_end[i])) {
    return(NA)
  } else {
    return(bmap$bin_end[i])
  }
}

ldscore <- ldscore %>% filter(CHR==chr) %>%
  rowwise() %>%
  mutate(bin_end=match_snp_to_bin(BP,bmap=bmap)) %>%
  mutate(CHR=as.numeric(CHR))

# Merge ldscore and bmap by bin
merged <- ldscore %>% left_join(bmap,by=c("bin_end","CHR"))

merged %>% select(CHR,SNP,BP,B) %>% fwrite(sep="\t", file=output_path)

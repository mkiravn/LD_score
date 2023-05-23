### File to write out and plot genome-wide correlation of LD score and B score
 
library(tidyverse)
library(data.table)

## Reading in files:
# LD score file
ld <- "data/UKBB.ALL.ldscore/UKBB.EUR.rsid.l2.ldscore.gz"
ld <- fread(ld,header=TRUE)

# function to match a SNP to its B score
match_snp_to_bin <- function(snp_bp, bmap) {
  i <- which.max(bmap$bin_end > snp_bp)
  if (i == 0 || is.na(bmap$bin_start[i]) || is.na(bmap$bin_end[i])) {
    return(NA)
  } else {
    return(bmap$bin_end[i])
  }
}

# Create a function to process each bmap file and join to ld score
process_bmap <- function(file,ldf=ld) {
  # Read the file and rename the columns
  bmap <- fread(file, header = TRUE) 
  colnames(bmap) <- c("B","length")
  chr_num <- as.numeric(gsub("[^0-9]+", "", file))
  # Add columns for bin start, bin end, and chromosome
  bmap <- bmap %>% 
    mutate(bin_end = cumsum(length), 
           bin_start = bin_end - length, 
           chr = chr_num)
  # join to ld score 
  ldb <- ldf %>%
    filter(CHR==chr_num) %>%
    rowwise() %>%
    mutate(bin_end=match_snp_to_bin(BP,bmap=bmap)) 

  ldb <- ldb %>% 
    rename("chr"="CHR","pos"="BP") %>%
    mutate(chr=as.numeric(chr),pos=as.numeric(pos)) %>%
    left_join(bmap,by=c("bin_end","chr")) %>%
    mutate(B=as.numeric(B)) %>%
    filter(!is.na(B)) %>%
    select(chr,pos,SNP,L2,B) 
    
  print(paste("Completed reading and joining file for chr",chr_num))
  return(ldb)
}


# Read all bmap files into a list and stack them
bmap_list <- lapply(list.files("data/bscores/CADD_bestfit", pattern = "\\.bmap\\.txt$", full.names = TRUE), process_bmap)

ld_bmap <- rbindlist(bmap_list)

# write out the file
write.table(ld_bmap,file="data/bscores/LD_Bscore.tsv",sep = "\t")



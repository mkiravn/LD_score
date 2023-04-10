library(tidyverse)
library(data.table)

bmap <- fread("../data/bscores/CADD_bestfit/chr1.bmap.txt")
old_bmap <- fread("../data/bscores/bkgd/chr1.bkgd")

colnames(bmap) <- c("B","length")
colnames(old_bmap) <- c("B","length")

bmap <- bmap %>% mutate(bin_end=cumsum(length),bin_start=bin_end-length,B=B/1000,CHR=1)
### PROBLEM: THE OLD ONES ARE ON HG18!!!! AND THE NEW ON HG19!!!
old_bmap <- old_bmap %>% mutate(bin_end=cumsum(length),bin_start=bin_end-length,CHR=1)

ldscore <- fread("../data/UKBB.ALL.ldscore/UKBB.EUR.rsid.l2.ldscore.gz")

sds <- fread("../data/sds/SDS_UK10K_n3195_release_Sep_19_2016.tab")

# Define function to match SNP BP to bin in bmap
match_snp_to_bin <- function(snp_bp, bmap) {
  i <- which.max(bmap$bin_end > snp_bp)
  if (i == 0 || is.na(bmap$bin_start[i]) || is.na(bmap$bin_end[i])) {
    return(NA)
  } else {
    return(bmap$bin_end[i])
  }
}

ldscore <- ldscore %>% filter(CHR==1) %>%
  rowwise() %>%
  mutate(bin_end=match_snp_to_bin(BP,bmap=bmap),
         bin_end_old=match_snp_to_bin(BP,bmap=old_bmap))

sds <- sds %>% filter(CHR==1) %>% rename("BP"="POS")

# Merge ldscore and bmap by bin
merged <- ldscore %>% left_join(bmap,by=c("bin_end","CHR")) %>%
#  left_join(old_bmap %>% rename("old_B"="B","bin_end_old"="bin_end"),by=c("bin_end_old","CHR")) %>%
  inner_join(sds,by=c("CHR","BP"))

# plot
merged %>% 
  ggplot(aes(y=L2,x=B)) +
  geom_point(alpha=0.1,col="grey") +
  geom_smooth(col="orchid") + 
  scale_y_log10() +
  theme_bw()

merged %>% 
  ggplot(aes(y=L2,x=old_B)) +
  geom_point(alpha=0.1,col="grey") +
  geom_smooth(col="orchid") + 
  scale_y_log10() +
  theme_bw()


merged %>% 
  ggplot(aes(y=B,x=old_B)) +
  geom_point(alpha=0.1,col="grey") +
  geom_smooth(col="orchid") + 
  theme_bw()

lm(L2 ~ B,data=merged) %>% summary()
lm(L2 ~ old_B,data=merged) %>% summary()

merged %>%
  group_by(bin_start,bin_end) %>%
  summarise(avg_L2=mean(L2),B=B) %>%
  ggplot(aes(x=B,y=avg_L2)) +
  geom_point(alpha=0.1,col="grey") +
  geom_smooth(col="orchid") + 
  scale_y_log10() +
  theme_bw()

# plot
merged %>% 
  ggplot(aes(y=SDS,x=B)) +
  geom_point(alpha=0.1,col="grey") +
  geom_smooth(col="orchid") + 
  theme_bw()

merged %>% 
  ggplot(aes(y=abs(SDS),x=B)) +
  geom_point(alpha=0.1,col="grey") +
  geom_smooth(col="orchid") + 
  theme_bw()

merged %>%
  group_by(bin_start,bin_end) %>%
  summarise(avg_SDS=mean(SDS),B=B) %>%
  ggplot(aes(x=B,y=avg_SDS)) +
  geom_point(alpha=0.1,col="grey") +
  geom_smooth(col="orchid") + 
  theme_bw()

# plot
merged %>% 
  ggplot(aes(x=SDS,y=B)) +
  geom_point(alpha=0.1,col="grey") +
  geom_smooth(col="orchid") + 
  theme_bw()

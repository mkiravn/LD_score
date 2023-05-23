library(data.table)
library(tidyverse)
sds <- "data/sds/30050.info0.allchroms.pan.tSDS.tsv"
sumstats <- "data/GWAS_summaries/processed/30050.info0.allchroms.pan.sumstats.tsv"
ld <- "data/UKBB.ALL.ldscore/UKBB.EUR.rsid.l2.ldscore.gz"
#bmap <- "data/bscores/CADD_bestfit/chr2.bmap.txt"

# first deal with ld
ld <- fread(ld,header=TRUE)
# now the others
tsds <- fread(sds,header=TRUE)
sumstats <- fread(sumstats,header=TRUE)
# Create a function to process each bmap file
process_bmap <- function(file) {
  # Read the file and rename the columns
  bmap <- fread(file, header = TRUE) 
  colnames(bmap) <- c("B","length")
  chr_num <- as.numeric(gsub("[^0-9]+", "", file))
  # Add columns for bin start, bin end, and chromosome
  bmap <- bmap %>% 
    mutate(bin_end = cumsum(length), 
           bin_start = bin_end - length, 
           chr = chr_num)
  print(paste("Reading bfile for chr",chr_num))
  return(bmap)
}

# Read all bmap files into a list and stack them
bmap_list <- lapply(list.files("data/bscores/CADD_bestfit", pattern = "\\.bmap\\.txt$", full.names = TRUE)[1:2], process_bmap)
bmap <- rbindlist(bmap_list)

match_snp_to_bin <- function(snp_bp, bmap) {
  i <- which.max(bmap$bin_end > snp_bp)
  if (i == 0 || is.na(bmap$bin_start[i]) || is.na(bmap$bin_end[i])) {
    return(NA)
  } else {
    return(bmap$bin_end[i])
  }
}
#print(head(bmap))
#print(head(ld))
#print(head(tsds))
#print(head(sumstats))
#print("processing LD")
ld <- ld %>%
  rowwise() %>%
  mutate(bin_end=match_snp_to_bin(BP,bmap=bmap))
print("joining")
joined <- inner_join(sumstats,tsds,by=c("chr","pos","SNP","A1","A2"))
#print(head(joined))
joined <- joined %>%
  mutate(cor=as.numeric(beta_EUR) * as.numeric(Z)) 
print(head(joined))
joined_all <- ld %>% 
  rename("chr"="CHR","pos"="BP") %>%
  mutate(chr=as.numeric(chr),pos=as.numeric(pos)) %>%
  left_join(bmap,by=c("bin_end","chr")) %>%
  right_join(joined,by=c("chr","pos","SNP")) %>% 
  mutate(B=as.numeric(B)) %>%
  filter(!is.na(B))
print(head(joined_all))

Bquantiles <- quantile(joined_all$B,na.rm=TRUE,probs = seq(0,1,0.1))
print(Bquantiles)

joined_all <- joined_all  %>%
  group_by(B) %>%
  mutate(Bquant=cut(B,breaks = c(0,Bquantiles,Inf))) %>%
  ungroup() 
#print(head(joined_all))
lm(cor~L2 + B,data=joined_all) %>% summary() %>% print()

joined_all %>%
  ggplot(aes(x=L2,y=cor)) +
  geom_point(aes(size=-log(pval),alpha=-log(pval)),col="tomato",shape=1) +
  geom_smooth(method="lm",col="black") +
  facet_wrap(~Bquant) +
  theme_bw() -> p1
  
lm(cor~L2 ,data=joined_all) %>% summary() %>% print()

joined_all %>% 
  mutate(rg=cor*L2) %>%
  group_by(Bquant) %>%
  summarise(mnrg=mean(rg),
            sdrg=sd(rg),
            mncor=mean(cor),
            sdcor=sd(cor),
    n=n()) %>%
  ggplot(aes(x=Bquant)) +
  geom_errorbar(aes(ymin=mnrg-sdrg*0.5,ymax=mnrg+sdrg*0.5)) +
  geom_point(aes(y=mnrg),size=3) +
  geom_line(aes(y=mnrg,group=1)) +
  geom_label(aes(y=1.5,x=Bquant,label=n),col="tomato4") +
  labs(y="mean  L2 x (tSDS x beta)",x="B score decile",
       title="All chromosomes") +
  theme_bw() -> p2

joined_all %>%
  group_by(Bquant) %>%
  summarise(zz=cor(beta_EUR,Z),
    rg=cor(cor,L2)) %>%
  ggplot(aes(x=Bquant,y=rg)) +
  geom_point(aes(y=zz),col="green") +
  geom_point() +
  theme_bw() -> p3


print("Saving plots")
dir.create("data/results/figs/")

ggsave(p1,filename=paste0("data/results/figs/","30050_p1.pdf"))
ggsave(p2,filename=paste0("data/results/figs/","30050_p2.pdf"))
ggsave(p3,filename=paste0("data/results/figs/","30050_p3.pdf"))

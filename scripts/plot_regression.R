library(data.table)
library(tidyverse)

tsds <- "data/sds/50_irnt.info0.chr2.tSDS.tsv"
sumstats <- "data/GWAS_summaries/processed/50_irnt.info0.chr2.sumstats.tsv"
ld <- "data/UKBB.ALL.ldscore/UKBB.EUR.rsid.l2.ldscore.gz"
bmap <- "data/bscores/CADD_bestfit/chr2.bmap.txt"

# first deal with ld
ld <- fread(ld,header=TRUE) %>% filter(CHR==2)
# now the others
tsds <- fread(tsds,header=TRUE)
sumstats <- fread(sumstats,header=TRUE)
bmap <- fread(bmap,header=TRUE)

colnames(bmap) <- c("B","length")

bmap <- bmap %>% mutate(bin_end=cumsum(length),bin_start=bin_end-length,chr=2)

match_snp_to_bin <- function(snp_bp, bmap) {
  i <- which.max(bmap$bin_end > snp_bp)
  if (i == 0 || is.na(bmap$bin_start[i]) || is.na(bmap$bin_end[i])) {
    return(NA)
  } else {
    return(bmap$bin_end[i])
  }
}

ld <- ld %>%
  rowwise() %>%
  mutate(bin_end=match_snp_to_bin(BP,bmap=bmap))

joined <- inner_join(sumstats,tsds,by=c("chr","pos","SNP","A1","A2")) %>%
  mutate(cor=beta * Z) 

joined_all <- ld %>% 
  rename("chr"="CHR","pos"="BP") %>%
  left_join(bmap,by=c("bin_end","chr")) %>%
  right_join(joined,by=c("chr","pos","SNP")) %>% 
  filter(!is.na(B))

Bquantiles <- quantile(joined_all$B,na.rm=TRUE,probs = seq(0,1,0.1))

joined_all <- joined_all  %>%
  group_by(B) %>%
  mutate(Bquant=cut(B,breaks = c(0,Bquantiles,Inf),labels = seq(0,1.1,0.1))) %>%
  ungroup() 

lm(cor~L2 + B,data=joined_all) %>% summary()

joined_all %>%
  ggplot(aes(x=L2,y=cor)) +
  geom_point(aes(size=-log(pval),alpha=-log(pval)),col="tomato",shape=1) +
  geom_smooth(method="lm",col="black") +
  #geom_label(aes(y=0,label=round(lm(cor~L2)$coefficients[2],4),x=0),col="tomato4") +
  facet_wrap(~Bquant) +
  theme_bw() 
  
lm(cor~L2 ,data=joined_all) %>% summary()

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
       title="Chromosome 2") +
  theme_bw()

# joined_all %>% 
#   mutate(rg=cor*L2) %>%
#   ggplot(aes(x=Bquant,y=rg)) +
#   geom_boxplot() +
#   theme_bw() +
#   labs(x="B score quantile")+
#   coord_flip() 

# 
# joined_all %>%
#   group_by(Bquant) %>%
#   summarise(rg=cor(cor,L2)) %>%
#   ggplot(aes(x=Bquant,y=rg)) +
#   geom_point() +
#   theme_bw() + 
#   theme(axis.text.x = element_text(angle = 90))

joined_all %>%
  group_by(Bquant) %>%
  summarise(zz=cor(beta,Z),
    rg=cor(cor,L2)) %>%
  ggplot(aes(x=Bquant,y=rg)) +
  geom_point(aes(y=zz),col="green") +
  geom_point() +
  theme_bw()

# joined_all %>%
#   mutate(L2=scale(L2)) %>%
#   mutate(zz=Z * beta,rg= zz * L2) %>%
#   ggplot(aes(x=Bquant,y=rg)) +
#   geom_boxplot(fill="tomato",col="tomato4",alpha=0.5,outlier.size = 0.1) +
#   #geom_boxplot(aes(y=zz),fill="green",col="darkgreen",alpha=0.5) +
#   theme_bw()

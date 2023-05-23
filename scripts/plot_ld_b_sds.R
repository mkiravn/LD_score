### LD and B score analysis with SDS

library(tidyverse)
library(data.table)

ld_b <- fread(file="data/bscores/LD_Bscore.tsv",header=T)

# calculate inverse B scores
ld_b <- ld_b %>% mutate(invB=1/B)

### now more analysis
# model of LD score as a function of B
ldb.model <- lm(L2~invB,data=ld_b) %>% summary() 
ldb.model %>% print()

ld_b %>% ggplot(aes(x=invB,y=L2)) +
  geom_point(size=0.1,col="lightgrey") +
  geom_density_2d(col="orchid") +
  geom_smooth(method="lm",col="black",se=F,lty=2) +
  theme_bw() +
  labs(caption=paste("r^2 =",round(ldb.model$r.squared,4)),
       x="1/B",
       y="L2") -> p1

## Separating quantiles of inverse B score
invBquantiles <- quantile(ld_b$invB,na.rm=TRUE,probs = seq(0,1,0.01))

ld_b.g <- ld_b %>% 
  group_by(invB) %>%
  mutate(invBquant=cut(invB,breaks = unique(invBquantiles))) %>%
  group_by(invBquant) %>%
  filter(!is.na(invBquant)) %>%
  summarise(mnL2=mean(L2),invBquant,mnInvB=mean(invB)) %>%
  ungroup() %>% 
  distinct() 

ld_b.g %>% ggplot(aes(x=mnInvB,y=mnL2)) +
  #geom_boxplot() +
  geom_point(col="orchid1") +
  geom_smooth(col="black",se=F,lty=2,size=0.5) +
  theme_bw() +
  labs(caption=paste("r^2 =",round(ldb.model$r.squared,4)),
       x="1/B (mean in quantile)",
       y="mean L2") -> p2

ggsave(
  plot = p1,
  path = "data/results/figs",
  device = "pdf",
  filename = "L2_B_1.pdf",
  width = 4,
  height = 3
)
ggsave(
  plot = p2,
  path = "data/results/figs",
  device = "pdf",
  filename = "L2_B_2.pdf",
  width = 4,
  height = 3
)


### Now bringing SDS into the picture
sds <- "data/sds/SDS_UK10K_n3195_release_Sep_19_2016.tab"
sds <- fread(file=sds,header=T)

ld_b_sds <- sds %>% 
  rename("chr"="CHR","pos"="POS","rsid"="ID") %>%
  inner_join(ld_b,by=c("chr","pos")) 

L2quantiles <- quantile(ld_b_sds$L2,na.rm=TRUE,probs = seq(0,1,0.01))

ld_b_sds.g <- ld_b_sds %>%
  rowwise() %>%
  mutate(L2quant=cut(L2,breaks = unique(L2quantiles))) %>%
  group_by(L2quant) %>%
  filter(!is.na(L2quant)) %>%
  summarise(mnL2=mean(L2),mnInvB=mean(invB),mnSDSchi2=mean(SDS^2)) %>%
  ungroup() %>% 
  distinct() 
  
ld_b_sds.g %>%   
  ggplot(aes(x=mnL2,y=mnSDSchi2)) +
  geom_point(col="orchid1") +
  geom_smooth(col="black",se=F,lty=2,size=0.5,method="lm") +
  theme_bw() +
  labs(x="L2 (mean in quantile)",
       y="SDS^2") -> p3

ld_b_sds.g %>%   
  ggplot(aes(x=mnInvB,y=mnSDSchi2)) +
  geom_point(col="orchid1") +
  geom_smooth(col="black",se=F,lty=2,size=0.5,method="lm") +
  theme_bw() +
  labs(x="1/B (mean in quantile)",
       y="SDS^2") -> p4

ggsave(
  plot = p3,
  path = "data/results/figs",
  device = "pdf",
  filename = "L2_B_SDS_1.pdf",
  width = 4,
  height = 3
)
ggsave(
  plot = p4,
  path = "data/results/figs",
  device = "pdf",
  filename = "L2_B_SDS_2.pdf",
  width = 4,
  height = 3
)



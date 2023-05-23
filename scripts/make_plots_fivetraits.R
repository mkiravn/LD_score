### A script to plot the relationship between tSDS * betahat and 
### LD score/B score for a variety of traits
library(tidyverse)
library(data.table)
library(ggpubr)

# specify the traits we're interested in
traits <- c("Standing Height",
            "Platelet crit",
            "Forced vital capacity (FVC)",
            "Trunk fat mass",
            "High light scatter reticulocyte percentage",
            "Pulse rate, automated reading",
            "Trunk predicted mass")

codes <- c(50,30090,3062,23128,30290,102,23130)

# read in the b scores
ld_b <- fread(file="data/bscores/LD_Bscore.tsv",header=T) %>%
		mutate(invB=1/B)

invBquantiles <- quantile(ld_b$invB,na.rm=TRUE,probs = seq(0,1,0.01))
L2quantiles <- quantile(ld_b$L2,na.rm=TRUE,probs = seq(0,1,0.01))

    ld_bg <- ld_b %>%
    group_by(chr) %>%
    mutate(wind = cut_interval(pos,10000)) %>%
    group_by(wind) %>%
    summarise(
      mnL2 = mean(L2),
      mnInvB = mean(invB),
     # mnrg = mean(rg),
     # mnbeta = mean(beta),
     # mntSDS=mean(Z),
      posn = min(pos),
      chr
	# lgpZ=max(-log10(pZ)),lgp=max(-log10(pval))
    ) %>%
    ungroup() %>%
    distinct() %>% arrange(chr) %>% mutate(strt=cumsum(as.numeric(posn)))


  ld_bg %>%
    ggplot(aes(x=strt)) +
      geom_point(aes(y=mnInvB,col=factor(chr)),size=0.5,alpha=0.5) +
      theme_bw() +
      #facet_wrap(~chr) +
      theme(axis.text.x = element_blank()) +
      geom_smooth(linewidth=0.5,se=F,aes(y=mnInvB,group=factor(chr)),col="black") +
      labs(x="position",y="1/B",col="chromosome") -> miamiB

  ld_bg %>%
    ggplot(aes(x=strt)) +
      geom_point(aes(y=mnL2,col=factor(chr)),size=0.5,alpha=0.5) +
      theme_bw() +
      #facet_wrap(~chr) +
      geom_smooth(linewidth=0.5,se=F,aes(y=mnL2,group=factor(chr)),col="black") +
      theme(axis.text.x = element_blank()) +
      labs(x="position",y="LD score",col="chromosome") -> miamiL

ggarrange(miamiB,miamiL,common.legend=T,nrow=2) -> miami

ggsave(plot=miami,"data/results/figs/corplots/miami.pdf",device="pdf",width=12,height=6)

make_dem_plots <- function(trait,code){
  
  filenames_sumstats <- file.path("data/GWAS_summaries/processed",
                                  paste0(code,".info0.allchroms.pan.sumstats.tsv"))
  filenames_tsds <- file.path("data/sds/",
                              paste0(code,".info0.allchroms.pan.tSDS.tsv"))
  
  
  sumstats <- fread(filenames_sumstats)
  tsds <- fread(filenames_tsds)
  rgldb <- inner_join(sumstats,tsds,
                      by=c("chr","pos","SNP","A1","A2")) %>%
    inner_join(ld_b,by=c("chr","pos")) %>% rename("beta"="beta_EUR")
   rgldb <- rgldb %>%  mutate(rg= beta * Z,
           pZ=pchisq(Z^2,1,lower.tail = F))

  rgldb.g <- rgldb %>%
    group_by(invB) %>%
    mutate(invBquant = cut(invB, breaks = unique(invBquantiles))) %>%
    group_by(invBquant) %>%
    filter(!is.na(invBquant)) %>%
    summarise(
      mnL2 = mean(L2),
      invBquant,
      mnInvB = mean(invB),
      mnrg = cor(beta,Z),
      mnbeta = mean(beta),
      mntSDS=mean(Z)
    ) %>%
    ungroup() %>%
    distinct() 
  
  rgldb.gl <- rgldb %>%
    group_by(L2) %>%
    mutate(L2quant = cut(L2, breaks = unique(L2quantiles))) %>%
    group_by(L2quant) %>%
    filter(!is.na(L2quant)) %>%
    summarise(
      mnL2 = mean(L2),
      L2quant,
      mnInvB = mean(invB),
      mnrg = cor(beta,Z),
      mnbeta = mean(beta),
      sdbeta = sd(beta),
      mntSDS=mean(Z),
      sdtSDS=sd(Z)
    ) %>%
    ungroup() %>%
    distinct()


  m <- lm(data=rgldb, rg ~ invB + L2)  %>% summary() 
  print(trait)
  print(m)
  m$coefficients[,c(1,2,4)] %>%
  as.data.frame() %>% rownames_to_column("v") %>%
  write_tsv(path=paste0("data/results/models/",code,"_model.txt"))
  
  rgldb.g %>%  
    ggplot(aes(x=mnInvB,y=mnrg))  +
    geom_hline(yintercept=0,col="black",lty=1) +
    geom_point(col="orchid1",size=0.5) +
    geom_smooth(method="lm",col="orchid4", lty=1,se=F,linewidth=0.5) +
    theme_bw() +
    geom_hline(yintercept=mean(rgldb$rg),col="orchid4",lty=2) +
    labs(x="1/B",y="Cor(tSDS,beta)",title=trait) -> p1
  
  rgldb.gl %>%  
    ggplot(aes(x=mnL2,y=mnrg))  +
    geom_hline(yintercept=0,col="black",lty=1) +
    geom_point(col="orchid1",size=0.5) +
    geom_smooth(method="lm",col="orchid4", lty=1,se=F,linewidth=0.5) +
    theme_bw() +
    geom_hline(yintercept=mean(rgldb.gl$mnrg),col="orchid4",lty=2) +
    labs(x="LD score",y="Cor(tSDS,beta)") -> p2
  
  rgldb.gl %>%  
    ggplot(aes(x=mnL2,y=mnbeta))  +
    geom_hline(yintercept=0,col="black",lty=1) +
    #geom_errorbar(aes(ymin=mnbeta-sdbeta,ymax=mnbeta+sdbeta))+
    geom_point(col="orchid1",size=0.5) +
    geom_smooth(method="lm",col="orchid4", lty=1,se=F,linewidth=0.5) +
    theme_bw() +
    labs(x="LD score",y="beta") +
    geom_hline(yintercept=mean(rgldb$beta),col="orchid4",lty=2) -> p3
  
  rgldb.gl %>%  
    ggplot(aes(x=mnL2,y=mntSDS))  +
    geom_hline(yintercept=0,col="black",lty=1) +
    #geom_errorbar(aes(ymin=mntSDS-sdtSDS,ymax=mntSDS+sdtSDS))+
    geom_point(col="orchid1",size=0.5) +
    geom_smooth(method="lm",col="orchid4", lty=1,se=F,linewidth=0.5) +
    theme_bw() +
    labs(x="LD score",y="tSDS") +
    geom_hline(yintercept=mean(tsds$Z),col="orchid4",lty=2) -> p4
    
    rgldb.gl %>%
    ggplot(aes(x=mnL2,y=sdtSDS^2))  +
    #geom_hline(yintercept=0,col="black",lty=1) +
    #geom_errorbar(aes(ymin=mntSDS-sdtSDS,ymax=mntSDS+sdtSDS))+
    geom_point(col="orchid1",size=0.5) +
    geom_smooth(method="lm",col="orchid4", lty=1,se=F,linewidth=0.5) +
    theme_bw() +
    labs(x="LD score",y="var(tSDS)") +
    geom_hline(yintercept=var(tsds$Z),col="orchid4",lty=2) -> p5
  
  rgldb.gl %>%
    ggplot(aes(x=mnL2,y=sdbeta^2))  +
    #geom_hline(yintercept=0,col="black",lty=1) +
    #geom_errorbar(aes(ymin=mntSDS-sdtSDS,ymax=mntSDS+sdtSDS))+
    geom_point(col="orchid1",size=0.5) +
    geom_smooth(method="lm",col="orchid4", lty=1,se=F,linewidth=0.5) +
    theme_bw() +
    labs(x="LD score",y="var(beta)") +
    geom_hline(yintercept=var(rgldb$beta),col="orchid4",lty=2) -> p6
  
  
  ggarrange(p1,p2,nrow=1) -> p12
  ggarrange(p3,p4,nrow=1) -> p34
  ggarrange(p5,p6,nrow=1) -> p56
  
  ggarrange(p34,p56,ncol=1) -> p3456
  
  ggarrange(p12,p3456,ncol=1,heights=c(1,2)) -> p123456
  ggsave(
    plot = p123456,
    height = 7,
    width = 5,
    filename = paste0(code, "_cors.pdf"),
    device = "pdf",
    path = "data/results/figs/corplots"
  )
 

}


for (i in c(1:length(traits))){
  trt <- traits[i]
  cd <- codes[i]
  make_dem_plots(trait = trt,code = cd)
}

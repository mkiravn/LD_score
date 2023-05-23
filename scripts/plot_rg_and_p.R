library(tidyverse)
library(stringr)
# Name pattern of the files
filename_pattern <- '*.info0.allchroms.*'



outpath <- "data/results/figs"
dir.create(outpath)


inpath <- "data/GWAS_summaries/processed"
# List all files in the directory
files <- list.files(inpath, pattern = filename_pattern, full.names = TRUE)
print(paste("Files:",files))

phenotypes <- read_tsv("data/simons_phenotypes_codes.tsv") %>% mutate(ncode=as.numeric(as.character(ncode))) 


for (filename in files){

	pheno_code <- str_extract(filename, "\\d+") %>% as.numeric()

#phenotypes <- read_tsv("data/simons_phenotypes_codes.tsv") %>%
#  mutate(ncode=as.numeric(as.character(ncode)))

	traitname <- phenotypes %>% filter(ncode==pheno_code) %>% select(Trait)

	tSDS <- read_tsv(file.path("data",
                           "sds",
                           paste0(pheno_code,
                                  ".info0.allchroms.pan.tSDS.tsv")))
	betas <- read_tsv(file.path(inpath,
                            paste0(pheno_code,
                                   ".info0.allchroms.pan.sumstats.tsv")))


	jn <- inner_join(tSDS,betas,by=c("SNP","chr","pos","A1","A2"))
	jn <- jn %>% mutate(pvalrank=rank(neglog10_pval_EUR),
                    pvalrankbin=cut(pvalrank,1000,ordered_result = T)) 

	jnsum <- jn %>% 
		group_by(pvalrankbin) %>%
		summarise(Z=mean(Z))

	cr <- cor(jn$Z,jn$pval) %>% round(4)

# ggplot(jn,aes(x=beta_EUR,y=Z)) +
#   geom_point(col="lightgrey",alpha=0.5) +
#   theme_bw()

# ggplot(jn,aes(x=pvalrank,y=Z)) +
#   geom_point(col="tomato") +
#   theme_bw() +
#   labs(x="p-value rank",y= "tSDS") +
#   geom_smooth()

ggplot(jnsum,aes(x=pvalrankbin,y=Z)) +
  geom_point(col="tomato") +
  theme_classic() +
  labs(x="p-value rank",y= "tSDS",subtitle = paste("Corr:",cr),
       title=traitname) +
  geom_smooth() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) -> pl

pl %>% ggsave(path=outpath,
              filename=paste0(pheno_code,
                              "rgplot.pdf"),
              device = "pdf",width = 4,height = 4,
              units = "in")
}

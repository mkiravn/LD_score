## Re-polarises SDS scores (AA/DA) to match UKBB (REF/AlT)

args=commandArgs(TRUE)


suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

sds_file = args[1]
ukbb_file = args[2]
sds_outfile = args[3]
ukbb_outfile = args[4]

####################
## Functions #######
####################
####################
SDS <- fread(sds_file) # in this file we have CHR POS ID AA DA
sumstats <-fread(ukbb_file) # in this we have chr pos id ref alt
# we want to join on chr pos id and switch the sign of SDS if AA!=ref and DA!=alt
# we could use this file to write a convention for the SDS and sumstats file
joined <- SDS %>% 
  rename("chr"="CHR","pos"="POS","rsid"="ID") %>%
  inner_join(sumstats,by=c("chr","pos","rsid")) 

print(paste("Signs to be flipped due to AA/ref differing:",length(joined$AA==joined$ref)))

joined <- joined %>%
  rowwise() %>%
  mutate(tSDS = case_when( 
    AA==ref & DA==alt ~ sign(beta) * SDS,
    AA==alt & DA==ref ~ - sign(beta) * SDS),
    .default=NA) %>% 
  ungroup() %>%
  filter(!is.na(tSDS))%>%
  mutate(N=3195.0) 

joined %>%
  select(chr,pos,SNP=rsid,A2,A1,Z=tSDS,N) %>%
  fwrite(sds_outfile, row.names = F, col.names = T, quote = F, sep = "\t")

joined %>%
  select(chr,pos,SNP=rsid,A2,A1,beta,pval,n_complete_samples) %>%
  fwrite(ukbb_outfile, row.names = F, col.names = T, quote = F, sep = "\t")
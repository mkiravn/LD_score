## Re-polarises SDS scores (AA/DA) to match UKBB (REF/AlT)

args=commandArgs(TRUE)


suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

sds_file = args[1]
ukbb_file = args[2]
outfile = args[3]

####################
## Functions #######
####################
####################
SDS <- fread(sds_file) # in this file we have CHR POS ID AA DA
sumstats <-fread(ukbb_file) # in this we have chr pos id ref alt
# we want to join on chr pos id and switch the sign of SDS if AA!=ref and DA!=alt

joined <- SDS %>% 
  rename("chr"="CHR","pos"="POS","rsid"="ID") %>%
  inner_join(sumstats,by=c("chr","pos","rsid"))

print(paste("Signs to be flipped due to AA/ref differing:",length(joined$AA==joined$ref)))

joined <- joined %>%
  rowwise() %>%
  mutate(tSDS = case_when( 
    AA==ref & DA==alt ~ sign(beta) * SDS,
    AA!=ref | DA!=alt ~ - sign(beta) * SDS)) %>% 
  ungroup() %>%
  mutate(N=3195.0) %>%
  select(chr,pos,SNP=rsid,A2=ref,A1=alt,Z=tSDS,N)

fwrite(joined, outfile, row.names = F, col.names = T, quote = F, sep = "\t")
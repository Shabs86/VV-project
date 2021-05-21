# AUTHOR - Shabbeer Hassan
# Date - 19/1/2021 
# Get independent loci from GWAS results

# Libraries ---------------------------------------------------------
setwd("F:/Dropbox")
getwd()
library("xlsx")
library("readxl")
library("ggplot2")
library("rsnps")
library("data.table")
library("ggrepel")
library("reshape2")
library("cowplot")
library("Rmisc")
library("dplyr")
library("stringr")
library("tidyr")
library("hexbin")
library("graphics")
library("readr")
#library("splines") # Used for the ns() function â (natural cubic splines)
library("tidyverse")
library("gapminder")
library("ggsci")
library("ggpubr")
library("extrafont")
library("sqldf")
#loadfonts(device = "win")
options(digits = 2)
options(scipen = 999)
library(annotables)

#install.packages(c("httr", "jsonlite"))
library(httr)
library(jsonlite)



# Input Dataset
VV_sig_only <- read_excel("Important_Documents/Doctoral_Work/Main Work/Paper 2 - Finngen_varicose veins_DVT/Variant lists etc/VV_sig_only.xlsx", 
                          sheet = "VV_Sig p-values only")

# Formatting dataset
#df = VV_sig_only %>% unite(Chr_Position, c(n, s), sep = " ", remove = FALSE)
#VV_summary_stats_lead_snps = mutate(VV_sig_only, concated_column = paste("chr", chrom, pos, sep = '_'))
VV_summary_stats_lead_snps = VV_sig_only %>% filter(VV_sig_only$pval <= 5e-8)
VV_summary_stats_lead_snps = mutate(VV_summary_stats_lead_snps, chr_variant = paste(chrom, pos,ref, alt, sep = ':'))

### Get LD correlations from API in form of a dataframe 

LD_server_API = "http://api.finngen.fi/api/ld?variant=blah&window=5000000&panel=sisu3&r2_thresh=0.1" #url for LD server
ld_api_all_chrom <- list() # Initiate list for final o/p

for(i in 1:nrow(VV_summary_stats_lead_snps)){
  ld_api <- vector(mode = "list", length = 0)
  ld_api = gsub("blah", VV_summary_stats_lead_snps$chr_variant[i], LD_server_API)
  VV_summary_stats_lead_snps$ld_api_url[i] <- ld_api
  #id <- VV_summary_stats_lead_snps$rsids[i]
  #list.name<- as.character(id)
  a <- GET(ld_api)
  b <- fromJSON(rawToChar(a$content))
  c <- data.frame(b$ld)
  ld_api_all_chrom <- rbind(ld_api_all_chrom, c)
}

# Add ld values into summary stat file
write.table(ld_api_all_chrom, file = "F:/Dropbox/Important_Documents/Doctoral_Work/Main Work/Paper 2 - Finngen_varicose veins_DVT/Variant lists etc/ld_api_all_chrom.txt",
            sep = "\t", row.names = FALSE, quote=F)
write.xlsx(ld_api_all_chrom, file = "F:/Dropbox/Important_Documents/Doctoral_Work/Main Work/Paper 2 - Finngen_varicose veins_DVT/Variant lists etc/ld_api_all_chrom.xlsx", sheetName = "Sheet_1",
            row.names = FALSE, col.names = TRUE)


# Write table into a text file
write.table(VV_summary_stats_lead_snps, file = "F:/Dropbox/Important_Documents/Doctoral_Work/Main Work/Paper 2 - Finngen_varicose veins_DVT/Variant lists etc/VV_summary_stats_lead_snps.txt",
            sep = "\t", row.names = FALSE, quote=F)



## Get pairwise LD values for snps in summary stat files
snps_with_ld_lead <- VV_summary_stats_lead_snps[VV_summary_stats_lead_snps$chr_variant %in% ld_api_all_chrom$variation2, ]
colnames(ld_api_all_chrom) <- c("d_prime", "r2", "variation1", "chr_variant")
LD_lead_snps <- inner_join(VV_summary_stats_lead_snps, ld_api_all_chrom, by = "chr_variant")
#LD_lead_snps <- merge(VV_summary_stats_lead_snps, ld_api_all_chrom, by = "chr_variant")
LD_lead_snps <- as.data.frame(LD_lead_snps)

# Write LD of lead snps into a text file
write.table(LD_lead_snps, file = "F:/Dropbox/Important_Documents/Doctoral_Work/Main Work/Paper 2 - Finngen_varicose veins_DVT/Variant lists etc/LD_lead_snps.txt",
            sep = "\t", row.names = FALSE, quote=F)
write.xlsx(LD_lead_snps, file = "F:/Dropbox/Important_Documents/Doctoral_Work/Main Work/Paper 2 - Finngen_varicose veins_DVT/Variant lists etc/LD_lead_snps.xlsx", sheetName = "Sheet_1",
           row.names = FALSE, col.names = TRUE)

## Import LD matrix
LD_matrix <- read_excel("Important_Documents/Doctoral_Work/Main Work/Paper 2 - Finngen_varicose veins_DVT/Variant lists etc/VV_summary_stats_lead_snps.xlsx", 
                        sheet = "LD_matrix", na = "NA")

## Get final Lead snps
To_be_removed_snps <- c("rs190947854", "rs191238524", "rs708461", "rs7977742", "rs201955556", "rs550822856", "rs571025773", "rs941510354", "rs182510184")
VV_summary_stats_lead_snps_with_LD <- as.data.frame(VV_summary_stats_lead_snps[!VV_summary_stats_lead_snps$rsids %in% To_be_removed_snps,])

# Write final Lead snps into a xlsx file
write.xlsx(VV_summary_stats_lead_snps_with_LD, file = "F:/Dropbox/Important_Documents/Doctoral_Work/Main Work/Paper 2 - Finngen_varicose veins_DVT/Variant lists etc/VV_summary_stats_lead_snps_with_LD.xlsx", sheetName = "Sheet_1",
           row.names = FALSE, col.names = TRUE)
write.table(VV_summary_stats_lead_snps_with_LD, file = "F:/Dropbox/Important_Documents/Doctoral_Work/Main Work/Paper 2 - Finngen_varicose veins_DVT/Variant lists etc/VV_summary_stats_lead_snps_with_LD.txt",
            sep = "\t", row.names = FALSE, quote=F)



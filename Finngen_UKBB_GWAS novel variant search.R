# AUTHOR - Shabbeer Hassan
# Date - 21/10/2020 
# Match variants across Finngen and UKBB 

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
library("splines") # Used for the ns() function â (natural cubic splines)
library("tidyverse")
library("gapminder")
library("ggsci")
library("ggpubr")
library("extrafont")
loadfonts(device = "win")
options(digits = 2)
options(scipen = 999)

# Input Arguments ---------------------------------------------------------

# Finngen_GWAS
Varicose_Veins_Autoreporting_Finngen_Gwas <- read_excel("Dropbox/Important_Documents/Doctoral_Work/Main Work/Paper 2 - Finngen_varicose veins_DVT/Variant lists etc/Varicose Veins_Autoreporting_Finngen_Gwas.xlsx", 
                                                        sheet = "Varicose Veins_Autoreporting_Fi")
View(Varicose_Veins_Autoreporting_Finngen_Gwas)
Below_cutoff <- data.frame(subset(Varicose_Veins_Autoreporting_Finngen_Gwas, pval< 5e-08))

loci <- unique(Below_cutoff$locus_id)
vars_list <- unique(Below_cutoff$X.variant)
credible_list <- unique(Below_cutoff$cs_id)

Below_cutoff %>% group_by(cs_id) %>% summarize(count=n()) 
Below_cutoff %>% group_by(X.variant) %>% summarize(count=n()) 
Below_cutoff %>% group_by(locus_id) %>% summarize(count=n()) 


# Merge coloc and finngen lists
coloc <- read_excel("Dropbox/Important_Documents/Doctoral_Work/Main Work/Paper 2 - Finngen_varicose veins_DVT/Variant lists etc/GTEx_coloc_results (1).xlsx")
View(coloc)

Final_list_Finngen <- read_excel("Dropbox/Important_Documents/Doctoral_Work/Main Work/Paper 2 - Finngen_varicose veins_DVT/Variant lists etc/Final list_Finngen.xlsx", 
                                 sheet = "Novel hits")
View(Final_list_Finngen)

All_sig_vars <- read_excel("Dropbox/Important_Documents/Doctoral_Work/Main Work/Paper 2 - Finngen_varicose veins_DVT/Variant lists etc/Final list_Finngen.xlsx", 
                           sheet = "All signif variants")
View(All_sig_vars)

Merged_all_sig_vars <- inner_join(All_sig_vars, coloc, 
                    by = c("Pos"))
Merged_novel_vars <- inner_join(Final_list_Finngen, coloc, 
                                  by = c("Pos"))

uniqie_loci_coloc <- unique(coloc$Locus)

# Output Files
write.xlsx(Merged_all_sig_vars, "Dropbox/Important_Documents/Doctoral_Work/Main Work/Paper 2 - Finngen_varicose veins_DVT/Variant lists etc/Merged_all_sig_vars.xlsx", sheetName = "Merged_all_sig_vars", 
  col.names = TRUE, row.names = TRUE, append = FALSE)

write.xlsx(Merged_novel_vars, "Dropbox/Important_Documents/Doctoral_Work/Main Work/Paper 2 - Finngen_varicose veins_DVT/Variant lists etc/Merged_all_sig_vars.xlsx", sheetName = "Merged_novel_vars", 
           col.names = TRUE, row.names = TRUE, append = TRUE)

## CONQUER Installation
install.packages("devtools")
install.packages("conquer.db_0.1.2.tar.gz", type="source", repos=NULL)
devtools::install_github("roderickslieker/CONQUER.DB")
devtools::install_github("roderickslieker/CONQUER.d3")
devtools::install_github("roderickslieker/CONQUER")


## CONQUER example

library(CONQUER)
DIR = "F:/Dropbox/Important_Documents/Doctoral_Work/Main Work/Paper 2 - Finngen_varicose veins_DVT"
library(readxl)
Final_list_Finngen <- read_excel("Important_Documents/Doctoral_Work/Main Work/Paper 2 - Finngen_varicose veins_DVT/Variant lists etc/Final list_Finngen.xlsx", 
                                 sheet = "Novel hits")

#snps <- Final_list_Finngen$RSID
snps <- c("rs192335285","rs191238524","rs635634","rs139988391","rs201955556") 

snps<-c("rs977371848", "rs992108547","rs947073006", "rs600038", "rs651007", "rs579459", "rs495828")


CONQUER::summarize(variants = snps ,
                   directory=DIR,
                   multiAnalyze=TRUE,
                   token='f3ebfb085237',
                   population = "FIN",
                   tissues=tissues_list)
#,
#                   tissues=c("Pancreas","Muscle_Skeletal","Liver"))

visualize(directory = DIR, SNPs = snps)

tissues_list <- conquer.db::gtexTissuesV8






if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenomicRanges")


install.packages("rjson_0.2.20.tar.gz", type="source", repos=NULL)


install.packages("enrichR_2.1.tar.gz", type="source", repos=NULL)


install.packages("shinyjs_2.0.0.tar.gz", type="source", repos=NULL)



# PRIMO analysis
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("limma")

install.packages("MASS")
install.packages("matrixStats")
install.packages("nnls")
install.packages("R.methodsS3")
install.packages("lcmix",repos="http://r-forge.r-project.org")

devtools::install_github("kjgleason/Primo")
library("Primo")


















############################# JUNK ##################################



# UKBB_GWAS
UKBB_gwas <- read_excel("Important_Documents/Doctoral_Work/Main Work/Paper 2 - Finngen_varicose veins_DVT/Variant lists etc/GWAS results from UKBB+23&me_2020 paper.xlsx")
View(UKBB_gwas)

# Replicated UKBB
UKBB_replicated <- read_excel("Important_Documents/Doctoral_Work/Main Work/Paper 2 - Finngen_varicose veins_DVT/Variant lists etc/Finngen_UKBB comp.xlsx", sheet = "Replicated_UKBB_23&Me")

# Unreplicated UKBB
UKBB_replicated <- read_excel("Important_Documents/Doctoral_Work/Main Work/Paper 2 - Finngen_varicose veins_DVT/Variant lists etc/Finngen_UKBB comp.xlsx", sheet = "Unreplicated_UKBB_23&Me")

# UKBB published
UKBB_published <- rbind(UKBB_replicated, UKBB_replicated)
colnames(UKBB_published)[which(names(UKBB_published) == "Position")] <- "pos"
colnames(UKBB_published)[which(names(UKBB_published) == "Chr")] <- "Chrom"

# Finngen_GWAS
Finngen <- read_excel("Important_Documents/Doctoral_Work/Main Work/Paper 2 - Finngen_varicose veins_DVT/Variant lists etc/Varicose Veins_Autoreporting_Finngen_Gwas.xlsx")
View(Finngen)

# Comparisons
UKBB_Finngen_common_3 <- Finngen %>% inner_join(UKBB_published, by=c("pos"))
a<-merge(UKBB_gwas, Finngen, by=c("pos"))



library(CONQUER)
snps <- c("rs11642430")

CONQUER::summarize(variants = snps,
                   directory=DIR,
                   multiAnalyze=TRUE,
                   token='f3ebfb085237',
                   tissues=c("Pancreas","Muscle_Skeletal","Liver"))
visualize(directory = DIR, SNPs = snps)

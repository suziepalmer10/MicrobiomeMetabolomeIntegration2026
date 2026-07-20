#dataset: Erawijantari
#purpose: data preparation
library("data.table")
library(dplyr)
library(tidyverse)
library(tidyr)
library(limma)
library(openxlsx)
library(writexl)

setwd("/home2/s180020/Desktop/IntegratedLearner_Revisions_02242026/A_Data_Processing/microbiome-metabolome-curated-data/data/processed_data/FRANZOSA_IBD_2019")


##load in datasets
#load in genera relative abundance
genera <- fread("genera.tsv", header = TRUE)
#load in species relative abundance
species <- fread("species.tsv", header=TRUE)
#load in metabolites
metabolites <- fread("mtb.tsv", header=TRUE) 
#metabolite annotations
mapped_metabolites <- fread("mtb.map.tsv", header=TRUE) 
#metadata
metadata <- fread("metadata.tsv", header=TRUE)

#SECTION 1: COMBINE, SCALE AND FILTER DATA. 
metadata_cols <- c("Sample", "Study.Group", 'Fecal.Calprotectin')

#check number of each condition
#make sure binary
counts <- table(metadata$Study.Group)
print(counts)

#add m__ before each metabolite. 
#this will allow us to parse and automate the data downstream
metabolites <- metabolites %>%
  rename_with(~ paste0("m__", .), -Sample)

#merge species and genus by Sample
mss <- merge(genera, species, by = "Sample")

#add t__ before all taxa to make easier to parse
filtered_mss<- mss %>%
  rename_with(~ paste0("t__", .), -Sample)


metadata_to_merge <- metadata[, ..metadata_cols]


#merge metadata and metabolites
metabolites_with_labels <- merge(metadata_to_merge, metabolites, by = "Sample")
processed_data <- merge(metabolites_with_labels, filtered_mss, by = "Sample")

#keep only control and CD
control_vs_cd_data <- subset(processed_data, Study.Group %in% c("Control", "CD"))
control_vs_uc_data <- subset(processed_data, Study.Group %in% c("Control", "UC"))
#convert CD and UC into IBD
control_vs_ibd <- processed_data
control_vs_ibd$Study.Group[control_vs_ibd$Study.Group %in% c("CD", "UC")] <- "IBD"


write.csv(control_vs_cd_data, file = "Franzosa_CD_data.csv", row.names = FALSE) 
write.csv(control_vs_uc_data, file = "Franzsosa_UC_data.csv", row.names = FALSE)
write.csv(control_vs_ibd, file = "Franzosa_IBD_data.csv", row.names = FALSE) 


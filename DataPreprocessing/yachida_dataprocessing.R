#dataset: Erawijantari
#purpose: data preparation
library("data.table")
library(dplyr)
library(tidyverse)
library(tidyr)
library(limma)
library(openxlsx)
library(writexl)

setwd("/home2/s180020/Desktop/IntegratedLearner_Revisions_02242026/A_Data_Processing/microbiome-metabolome-curated-data/data/processed_data/YACHIDA_CRC_2019")


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
metadata_cols <- c("Sample", "Study.Group")

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

#keep only control and stage I_II
control_vs_stage1_2 <- subset(processed_data, Study.Group %in% c("Healthy", "Stage_I_II"))
control_vs_stage3_4 <- subset(processed_data, Study.Group %in% c("Healthy", "Stage_III_IV"))
#control vs MP Stage 0
control_vs_MP_stage0<- subset(processed_data, Study.Group %in% c("Healthy", "Stage_0", 'MP'))
control_vs_MP_stage0$Study.Group[control_vs_MP_stage0$Study.Group %in% c("Stage_0", "MP")] <- "MP_stage0"
#control vs all CRC
control_vs_CRC<- subset(processed_data, Study.Group %in% c("Healthy", "Stage_0", 'MP','Stage_I_II', 'Stage_III_IV'))
control_vs_CRC$Study.Group[control_vs_CRC$Study.Group %in% c("Stage_0", 'MP','Stage_I_II', 'Stage_III_IV')] <- "CRC"

counts <- table(control_vs_stage1_2$Study.Group)
print(counts)
counts <- table(control_vs_stage3_4$Study.Group)
print(counts)
counts <- table(control_vs_MP_stage0$Study.Group)
print(counts)
counts <- table(control_vs_CRC$Study.Group)
print(counts)




write.csv(control_vs_stage1_2, file = "Yachida_control_vs_stage1_2.csv", row.names = FALSE) 
write.csv(control_vs_stage3_4, file = "Yachida_control_vs_stage3_4.csv", row.names = FALSE) 
write.csv(control_vs_MP_stage0, file = "Yachida_control_vs_MP_stage0.csv", row.names = FALSE) 
write.csv(control_vs_CRC, file = "Yachida_control_vs_CRC.csv", row.names = FALSE) 


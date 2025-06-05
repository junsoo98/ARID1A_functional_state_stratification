library("KEGGREST")
library("progress")
library("tidyverse") 
library("dplyr")
library("TCGAbiolinks") 
library("stringr") 
library("TCGAutils") 
library('tibble') 
library('caret') 
library('datawizard') 
library(factoextra)
library(FactoMineR)

######################################
# KEGG pathway database
# - hsa +
# * Transcription
# 03020 RNA polymerase *
# 03022 Basal transcription factors *
# 03040 Spliceosome
# * Translation
# 03010 Ribosome
# 00970 Aminoacyl-tRNA biosynthesis
# 03013 Nucleocytoplasmic transport
# 03015 mRNA surveillance pathway
# 03008 Ribosome biogenesis in eukaryotes
# * Degradation
# 03060 Protein export
# 04141 Protein processing in endoplasmic reticulum
# 04130 SNARE interactions in vesicular transport
# 04120 Ubiquitin mediated proteolysis
# 04122 Sulfur relay system
# 03050 Proteasome
# 03018 RNA degradation
#######################################

RNA_polymerase <-  keggGet("hsa03020")[[1]]$GENE
RNA_polymerase <-  RNA_polymerase[seq(0,length(RNA_polymerase),2)]
RNA_polymerase <- gsub("\\;.*","",RNA_polymerase)

Basal_transcription_factors <-  keggGet("hsa03022")[[1]]$GENE
Basal_transcription_factors <-  Basal_transcription_factors[seq(0,length(Basal_transcription_factors),2)]
Basal_transcription_factors <- gsub("\\;.*","",Basal_transcription_factors)

Spliceosome <-  keggGet("hsa03040")[[1]]$GENE
Spliceosome <-  Spliceosome[seq(0,length(Spliceosome),2)]
Spliceosome <- gsub("\\;.*","",Spliceosome)

Ribosome <-  keggGet("hsa03010")[[1]]$GENE
Ribosome <-  Ribosome[seq(0,length(Ribosome),2)]
Ribosome <- gsub("\\;.*","",Ribosome)

Aminoacyl_tRNA_biosynthesis <-  keggGet("hsa00970")[[1]]$GENE
Aminoacyl_tRNA_biosynthesis <-  Aminoacyl_tRNA_biosynthesis[seq(0,length(Aminoacyl_tRNA_biosynthesis),2)]
Aminoacyl_tRNA_biosynthesis <- gsub("\\;.*","",Aminoacyl_tRNA_biosynthesis)

Nucleocytoplasmic_transport <-  keggGet("hsa03013")[[1]]$GENE
Nucleocytoplasmic_transport <-  Nucleocytoplasmic_transport[seq(0,length(Nucleocytoplasmic_transport),2)]
Nucleocytoplasmic_transport <- gsub("\\;.*","",Nucleocytoplasmic_transport)

mRNA_surveillance_pathway <-  keggGet("hsa03015")[[1]]$GENE
mRNA_surveillance_pathway <-  mRNA_surveillance_pathway[seq(0,length(mRNA_surveillance_pathway),2)]
mRNA_surveillance_pathway <- gsub("\\;.*","",mRNA_surveillance_pathway)

Ribosome_biogenesis_in_eukaryotes <-  keggGet("hsa03008")[[1]]$GENE
Ribosome_biogenesis_in_eukaryotes <-  Ribosome_biogenesis_in_eukaryotes[seq(0,length(Ribosome_biogenesis_in_eukaryotes),2)]
Ribosome_biogenesis_in_eukaryotes <- gsub("\\;.*","",Ribosome_biogenesis_in_eukaryotes)

Protein_export <-  keggGet("hsa03060")[[1]]$GENE
Protein_export <-  Protein_export[seq(0,length(Protein_export),2)]
Protein_export <- gsub("\\;.*","",Protein_export)

Protein_processing_in_endoplasmic_reticulum <-  keggGet("hsa04141")[[1]]$GENE
Protein_processing_in_endoplasmic_reticulum <-  Protein_processing_in_endoplasmic_reticulum[seq(0,length(Protein_processing_in_endoplasmic_reticulum),2)]
Protein_processing_in_endoplasmic_reticulum <- gsub("\\;.*","",Protein_processing_in_endoplasmic_reticulum)

SNARE_interactions <-  keggGet("hsa04130")[[1]]$GENE
SNARE_interactions <-  SNARE_interactions[seq(0,length(SNARE_interactions),2)]
SNARE_interactions <- gsub("\\;.*","",SNARE_interactions)

Ubiquitin <-  keggGet("hsa04120")[[1]]$GENE
Ubiquitin <-  Ubiquitin[seq(0,length(Ubiquitin),2)]
Ubiquitin <- gsub("\\;.*","",Ubiquitin)

Sulfur_relay_system <-  keggGet("hsa04122")[[1]]$GENE
Sulfur_relay_system <-  Sulfur_relay_system[seq(0,length(Sulfur_relay_system),2)]
Sulfur_relay_system <- gsub("\\;.*","",Sulfur_relay_system)

Proteasome <-  keggGet("hsa03050")[[1]]$GENE
Proteasome <-  Proteasome[seq(0,length(Proteasome),2)]
Proteasome <- gsub("\\;.*","",Proteasome)

RNA_degradation <-  keggGet("hsa03018")[[1]]$GENE
RNA_degradation <-  RNA_degradation[seq(0,length(RNA_degradation),2)]
RNA_degradation <- gsub("\\;.*","",RNA_degradation)

Pathway_list <- c("RNA_polymerase", "Basal_transcription_factors", "Spliceosome", "Ribosome","Aminoacyl_tRNA_biosynthesis",
                  "Nucleocytoplasmic_transport","mRNA_surveillance_pathway","Ribosome_biogenesis_in_eukaryotes",
                  "Protein_export","Protein_processing_in_endoplasmic_reticulum","SNARE_interactions",
                  "Ubiquitin", "Sulfur_relay_system", "Proteasome", "RNA_degradation")

RNA_sample_sheet_path <- paste0(data_dir, "/RNA_sample_sheet.tsv")
RNA_sample_sheet <- read.delim(RNA_sample_sheet_path, sep = "\t")
RNA_sample_sheet <- RNA_sample_sheet %>% 
  filter(Sample.Type == "Primary Tumor") %>% 
  filter(str_detect(Sample.ID, "-01A"))


#####################
## mRNA processing ##
#####################
data_dir <- "dataset"
RNA_dir <- paste0(data_dir,"/RNA_seq")
All_RNA_dir <- list.dirs(RNA_dir, recursive = FALSE)

#####
RNA_polymerase_frame <- data.frame(matrix(nrow = length(unique(RNA_sample_sheet[,"Case.ID"])), 
                               ncol = length(RNA_polymerase)))
rownames(RNA_polymerase_frame) <- unique(RNA_sample_sheet[,"Case.ID"])
colnames(RNA_polymerase_frame) <- RNA_polymerase
RNA_polymerase_frame[] <- NaN
pb <- progress_bar$new(total = length(All_RNA_dir)) #checking for loop progress
for (dir_name in All_RNA_dir) {
  pb$tick()
  rna_file_name <- list.files(dir_name, pattern = "*.tsv")
  if (!(rna_file_name %in% RNA_sample_sheet$File.Name)) {
    next
  }
  rna_file <- read.delim(paste0(dir_name,"/",rna_file_name),
                         sep = '\t',header = T,fill = T, skip = 1)
  for (component in RNA_polymerase) {
    if (!(component %in% rna_file$gene_name)) {
      next
    }
    if (is.nan(RNA_polymerase_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],component])) { ##
      RNA_polymerase_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],component] <- sum(rna_file[rna_file$gene_name == component,"unstranded"])
    }
    else {
      cur <- RNA_polymerase_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],component]
      RNA_polymerase_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],component] <- round(cur + sum(rna_file[rna_file$gene_name == component,"unstranded"])/2)
    }
  }
}
#remove non-data gene
RNA_polymerase_frame <- RNA_polymerase_frame[,!is.nan(colSums(RNA_polymerase_frame))]
#feature scaling
process <- preProcess(RNA_polymerase_frame, method = c("range"))
scaled_RNA_polymerase_frame <- predict(process, RNA_polymerase_frame)

#####
Basal_transcription_factors_frame <- data.frame(matrix(nrow = length(unique(RNA_sample_sheet[,"Case.ID"])), 
                                          ncol = length(Basal_transcription_factors)))
rownames(Basal_transcription_factors_frame) <- unique(RNA_sample_sheet[,"Case.ID"])
colnames(Basal_transcription_factors_frame) <- Basal_transcription_factors
Basal_transcription_factors_frame[] <- NaN
pb <- progress_bar$new(total = length(All_RNA_dir)) #checking for loop progress
for (dir_name in All_RNA_dir) {
  pb$tick()
  rna_file_name <- list.files(dir_name, pattern = "*.tsv")
  if (!(rna_file_name %in% RNA_sample_sheet$File.Name)) {
    next
  }
  rna_file <- read.delim(paste0(dir_name,"/",rna_file_name),
                         sep = '\t',header = T,fill = T, skip = 1)
  for (component in Basal_transcription_factors) {
    if (!(component %in% rna_file$gene_name)) {
      next
    }
    if (is.nan(Basal_transcription_factors_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],component])) { ##
      Basal_transcription_factors_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],component] <- sum(rna_file[rna_file$gene_name == component,"unstranded"])
    }
    else {
      cur <- Basal_transcription_factors_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],component]
      Basal_transcription_factors_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],component] <- round(cur + sum(rna_file[rna_file$gene_name == component,"unstranded"])/2)
    }
  }
}
#remove non-data gene
Basal_transcription_factors_frame <- Basal_transcription_factors_frame[,!is.nan(colSums(Basal_transcription_factors_frame))]
#feature scaling
process <- preProcess(Basal_transcription_factors_frame, method = c("range"))
scaled_Basal_transcription_factors_frame <- predict(process, Basal_transcription_factors_frame)

#####
Spliceosome_frame <- data.frame(matrix(nrow = length(unique(RNA_sample_sheet[,"Case.ID"])), 
                                          ncol = length(Spliceosome)))
rownames(Spliceosome_frame) <- unique(RNA_sample_sheet[,"Case.ID"])
colnames(Spliceosome_frame) <- Spliceosome
Spliceosome_frame[] <- NaN
pb <- progress_bar$new(total = length(All_RNA_dir)) #checking for loop progress
for (dir_name in All_RNA_dir) {
  pb$tick()
  rna_file_name <- list.files(dir_name, pattern = "*.tsv")
  if (!(rna_file_name %in% RNA_sample_sheet$File.Name)) {
    next
  }
  rna_file <- read.delim(paste0(dir_name,"/",rna_file_name),
                         sep = '\t',header = T,fill = T, skip = 1)
  for (component in Spliceosome) {
    if (!(component %in% rna_file$gene_name)) {
      next
    }
    if (is.nan(Spliceosome_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],component])) { ##
      Spliceosome_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],component] <- sum(rna_file[rna_file$gene_name == component,"unstranded"])
    }
    else {
      cur <- Spliceosome_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],component]
      Spliceosome_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],component] <- round(cur + sum(rna_file[rna_file$gene_name == component,"unstranded"])/2)
    }
  }
}
#remove non-data gene
Spliceosome_frame <- Spliceosome_frame[,!is.nan(colSums(Spliceosome_frame))]
#feature scaling
process <- preProcess(Spliceosome_frame, method = c("range"))
scaled_Spliceosome_frame <- predict(process, Spliceosome_frame)

#####
Ribosome_frame <- data.frame(matrix(nrow = length(unique(RNA_sample_sheet[,"Case.ID"])), 
                                          ncol = length(Ribosome)))
rownames(Ribosome_frame) <- unique(RNA_sample_sheet[,"Case.ID"])
colnames(Ribosome_frame) <- Ribosome
Ribosome_frame[] <- NaN
pb <- progress_bar$new(total = length(All_RNA_dir)) #checking for loop progress
for (dir_name in All_RNA_dir) {
  pb$tick()
  rna_file_name <- list.files(dir_name, pattern = "*.tsv")
  if (!(rna_file_name %in% RNA_sample_sheet$File.Name)) {
    next
  }
  rna_file <- read.delim(paste0(dir_name,"/",rna_file_name),
                         sep = '\t',header = T,fill = T, skip = 1)
  for (component in Ribosome) {
    if (!(component %in% rna_file$gene_name)) {
      next
    }
    if (is.nan(Ribosome_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],component])) { ##
      Ribosome_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],component] <- sum(rna_file[rna_file$gene_name == component,"unstranded"])
    }
    else {
      cur <- Ribosome_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],component]
      Ribosome_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],component] <- round(cur + sum(rna_file[rna_file$gene_name == component,"unstranded"])/2)
    }
  }
}
#remove non-data gene
Ribosome_frame <- Ribosome_frame[,!is.nan(colSums(Ribosome_frame))]
#feature scaling
process <- preProcess(Ribosome_frame, method = c("range"))
scaled_Ribosome_frame <- predict(process, Ribosome_frame)

#####
Aminoacyl_tRNA_biosynthesis_frame <- data.frame(matrix(nrow = length(unique(RNA_sample_sheet[,"Case.ID"])), 
                                          ncol = length(Aminoacyl_tRNA_biosynthesis)))
rownames(Aminoacyl_tRNA_biosynthesis_frame) <- unique(RNA_sample_sheet[,"Case.ID"])
colnames(Aminoacyl_tRNA_biosynthesis_frame) <- Aminoacyl_tRNA_biosynthesis
Aminoacyl_tRNA_biosynthesis_frame[] <- NaN
pb <- progress_bar$new(total = length(All_RNA_dir)) #checking for loop progress
for (dir_name in All_RNA_dir) {
  pb$tick()
  rna_file_name <- list.files(dir_name, pattern = "*.tsv")
  if (!(rna_file_name %in% RNA_sample_sheet$File.Name)) {
    next
  }
  rna_file <- read.delim(paste0(dir_name,"/",rna_file_name),
                         sep = '\t',header = T,fill = T, skip = 1)
  for (component in Aminoacyl_tRNA_biosynthesis) {
    if (!(component %in% rna_file$gene_name)) {
      next
    }
    if (is.nan(Aminoacyl_tRNA_biosynthesis_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],component])) { ##
      Aminoacyl_tRNA_biosynthesis_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],component] <- sum(rna_file[rna_file$gene_name == component,"unstranded"])
    }
    else {
      cur <- Aminoacyl_tRNA_biosynthesis_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],component]
      Aminoacyl_tRNA_biosynthesis_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],component] <- round(cur + sum(rna_file[rna_file$gene_name == component,"unstranded"])/2)
    }
  }
}
#remove non-data gene
Aminoacyl_tRNA_biosynthesis_frame <- Aminoacyl_tRNA_biosynthesis_frame[,!is.nan(colSums(Aminoacyl_tRNA_biosynthesis_frame))]
#feature scaling
process <- preProcess(Aminoacyl_tRNA_biosynthesis_frame, method = c("range"))
scaled_Aminoacyl_tRNA_biosynthesis_frame <- predict(process, Aminoacyl_tRNA_biosynthesis_frame)

#####
Nucleocytoplasmic_transport_frame <- data.frame(matrix(nrow = length(unique(RNA_sample_sheet[,"Case.ID"])), 
                                          ncol = length(Nucleocytoplasmic_transport)))
rownames(Nucleocytoplasmic_transport_frame) <- unique(RNA_sample_sheet[,"Case.ID"])
colnames(Nucleocytoplasmic_transport_frame) <- Nucleocytoplasmic_transport
Nucleocytoplasmic_transport_frame[] <- NaN
pb <- progress_bar$new(total = length(All_RNA_dir)) #checking for loop progress
for (dir_name in All_RNA_dir) {
  pb$tick()
  rna_file_name <- list.files(dir_name, pattern = "*.tsv")
  if (!(rna_file_name %in% RNA_sample_sheet$File.Name)) {
    next
  }
  rna_file <- read.delim(paste0(dir_name,"/",rna_file_name),
                         sep = '\t',header = T,fill = T, skip = 1)
  for (component in Nucleocytoplasmic_transport) {
    if (!(component %in% rna_file$gene_name)) {
      next
    }
    if (is.nan(Nucleocytoplasmic_transport_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],component])) { ##
      Nucleocytoplasmic_transport_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],component] <- sum(rna_file[rna_file$gene_name == component,"unstranded"])
    }
    else {
      cur <- Nucleocytoplasmic_transport_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],component]
      Nucleocytoplasmic_transport_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],component] <- round(cur + sum(rna_file[rna_file$gene_name == component,"unstranded"])/2)
    }
  }
}
#remove non-data gene
Nucleocytoplasmic_transport_frame <- Nucleocytoplasmic_transport_frame[,!is.nan(colSums(Nucleocytoplasmic_transport_frame))]
#feature scaling
process <- preProcess(Nucleocytoplasmic_transport_frame, method = c("range"))
scaled_Nucleocytoplasmic_transport_frame <- predict(process, Nucleocytoplasmic_transport_frame)

#####
mRNA_surveillance_pathway_frame <- data.frame(matrix(nrow = length(unique(RNA_sample_sheet[,"Case.ID"])), 
                                          ncol = length(mRNA_surveillance_pathway)))
rownames(mRNA_surveillance_pathway_frame) <- unique(RNA_sample_sheet[,"Case.ID"])
colnames(mRNA_surveillance_pathway_frame) <- mRNA_surveillance_pathway
mRNA_surveillance_pathway_frame[] <- NaN
pb <- progress_bar$new(total = length(All_RNA_dir)) #checking for loop progress
for (dir_name in All_RNA_dir) {
  pb$tick()
  rna_file_name <- list.files(dir_name, pattern = "*.tsv")
  if (!(rna_file_name %in% RNA_sample_sheet$File.Name)) {
    next
  }
  rna_file <- read.delim(paste0(dir_name,"/",rna_file_name),
                         sep = '\t',header = T,fill = T, skip = 1)
  for (component in mRNA_surveillance_pathway) {
    if (!(component %in% rna_file$gene_name)) {
      next
    }
    if (is.nan(mRNA_surveillance_pathway_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],component])) { ##
      mRNA_surveillance_pathway_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],component] <- sum(rna_file[rna_file$gene_name == component,"unstranded"])
    }
    else {
      cur <- mRNA_surveillance_pathway_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],component]
      mRNA_surveillance_pathway_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],component] <- round(cur + sum(rna_file[rna_file$gene_name == component,"unstranded"])/2)
    }
  }
}
#remove non-data gene
mRNA_surveillance_pathway_frame <- mRNA_surveillance_pathway_frame[,!is.nan(colSums(mRNA_surveillance_pathway_frame))]
#feature scaling
process <- preProcess(mRNA_surveillance_pathway_frame, method = c("range"))
scaled_mRNA_surveillance_pathway_frame <- predict(process, mRNA_surveillance_pathway_frame)

#####
Ribosome_biogenesis_in_eukaryotes_frame <- data.frame(matrix(nrow = length(unique(RNA_sample_sheet[,"Case.ID"])), 
                                          ncol = length(Ribosome_biogenesis_in_eukaryotes)))
rownames(Ribosome_biogenesis_in_eukaryotes_frame) <- unique(RNA_sample_sheet[,"Case.ID"])
colnames(Ribosome_biogenesis_in_eukaryotes_frame) <- Ribosome_biogenesis_in_eukaryotes
Ribosome_biogenesis_in_eukaryotes_frame[] <- NaN
pb <- progress_bar$new(total = length(All_RNA_dir)) #checking for loop progress
for (dir_name in All_RNA_dir) {
  pb$tick()
  rna_file_name <- list.files(dir_name, pattern = "*.tsv")
  if (!(rna_file_name %in% RNA_sample_sheet$File.Name)) {
    next
  }
  rna_file <- read.delim(paste0(dir_name,"/",rna_file_name),
                         sep = '\t',header = T,fill = T, skip = 1)
  for (component in Ribosome_biogenesis_in_eukaryotes) {
    if (!(component %in% rna_file$gene_name)) {
      next
    }
    if (is.nan(Ribosome_biogenesis_in_eukaryotes_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],component])) { ##
      Ribosome_biogenesis_in_eukaryotes_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],component] <- sum(rna_file[rna_file$gene_name == component,"unstranded"])
    }
    else {
      cur <- Ribosome_biogenesis_in_eukaryotes_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],component]
      Ribosome_biogenesis_in_eukaryotes_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],component] <- round(cur + sum(rna_file[rna_file$gene_name == component,"unstranded"])/2)
    }
  }
}
#remove non-data gene
Ribosome_biogenesis_in_eukaryotes_frame <- Ribosome_biogenesis_in_eukaryotes_frame[,!is.nan(colSums(Ribosome_biogenesis_in_eukaryotes_frame))]
#feature scaling
process <- preProcess(Ribosome_biogenesis_in_eukaryotes_frame, method = c("range"))
scaled_Ribosome_biogenesis_in_eukaryotes_frame <- predict(process, Ribosome_biogenesis_in_eukaryotes_frame)

#####
Protein_export_frame <- data.frame(matrix(nrow = length(unique(RNA_sample_sheet[,"Case.ID"])), 
                                          ncol = length(Protein_export)))
rownames(Protein_export_frame) <- unique(RNA_sample_sheet[,"Case.ID"])
colnames(Protein_export_frame) <- Protein_export
Protein_export_frame[] <- NaN
pb <- progress_bar$new(total = length(All_RNA_dir)) #checking for loop progress
for (dir_name in All_RNA_dir) {
  pb$tick()
  rna_file_name <- list.files(dir_name, pattern = "*.tsv")
  if (!(rna_file_name %in% RNA_sample_sheet$File.Name)) {
    next
  }
  rna_file <- read.delim(paste0(dir_name,"/",rna_file_name),
                         sep = '\t',header = T,fill = T, skip = 1)
  for (component in Protein_export) {
    if (!(component %in% rna_file$gene_name)) {
      next
    }
    if (is.nan(Protein_export_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],component])) { ##
      Protein_export_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],component] <- sum(rna_file[rna_file$gene_name == component,"unstranded"])
    }
    else {
      cur <- Protein_export_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],component]
      Protein_export_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],component] <- round(cur + sum(rna_file[rna_file$gene_name == component,"unstranded"])/2)
    }
  }
}
#remove non-data gene
Protein_export_frame <- Protein_export_frame[,!is.nan(colSums(Protein_export_frame))]
#feature scaling
process <- preProcess(Protein_export_frame, method = c("range"))
scaled_Protein_export_frame <- predict(process, Protein_export_frame)

#####
Protein_processing_in_endoplasmic_reticulum_frame <- data.frame(matrix(nrow = length(unique(RNA_sample_sheet[,"Case.ID"])), 
                                          ncol = length(Protein_processing_in_endoplasmic_reticulum)))
rownames(Protein_processing_in_endoplasmic_reticulum_frame) <- unique(RNA_sample_sheet[,"Case.ID"])
colnames(Protein_processing_in_endoplasmic_reticulum_frame) <- Protein_processing_in_endoplasmic_reticulum
Protein_processing_in_endoplasmic_reticulum_frame[] <- NaN
pb <- progress_bar$new(total = length(All_RNA_dir)) #checking for loop progress
for (dir_name in All_RNA_dir) {
  pb$tick()
  rna_file_name <- list.files(dir_name, pattern = "*.tsv")
  if (!(rna_file_name %in% RNA_sample_sheet$File.Name)) {
    next
  }
  rna_file <- read.delim(paste0(dir_name,"/",rna_file_name),
                         sep = '\t',header = T,fill = T, skip = 1)
  for (component in Protein_processing_in_endoplasmic_reticulum) {
    if (!(component %in% rna_file$gene_name)) {
      next
    }
    if (is.nan(Protein_processing_in_endoplasmic_reticulum_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],component])) { ##
      Protein_processing_in_endoplasmic_reticulum_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],component] <- sum(rna_file[rna_file$gene_name == component,"unstranded"])
    }
    else {
      cur <- Protein_processing_in_endoplasmic_reticulum_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],component]
      Protein_processing_in_endoplasmic_reticulum_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],component] <- round(cur + sum(rna_file[rna_file$gene_name == component,"unstranded"])/2)
    }
  }
}
#remove non-data gene
Protein_processing_in_endoplasmic_reticulum_frame <- Protein_processing_in_endoplasmic_reticulum_frame[,!is.nan(colSums(Protein_processing_in_endoplasmic_reticulum_frame))]
#feature scaling
process <- preProcess(Protein_processing_in_endoplasmic_reticulum_frame, method = c("range"))
scaled_Protein_processing_in_endoplasmic_reticulum_frame <- predict(process, Protein_processing_in_endoplasmic_reticulum_frame)

#####
SNARE_interactions_frame <- data.frame(matrix(nrow = length(unique(RNA_sample_sheet[,"Case.ID"])), 
                                          ncol = length(SNARE_interactions)))
rownames(SNARE_interactions_frame) <- unique(RNA_sample_sheet[,"Case.ID"])
colnames(SNARE_interactions_frame) <- SNARE_interactions
SNARE_interactions_frame[] <- NaN
pb <- progress_bar$new(total = length(All_RNA_dir)) #checking for loop progress
for (dir_name in All_RNA_dir) {
  pb$tick()
  rna_file_name <- list.files(dir_name, pattern = "*.tsv")
  if (!(rna_file_name %in% RNA_sample_sheet$File.Name)) {
    next
  }
  rna_file <- read.delim(paste0(dir_name,"/",rna_file_name),
                         sep = '\t',header = T,fill = T, skip = 1)
  for (component in SNARE_interactions) {
    if (!(component %in% rna_file$gene_name)) {
      next
    }
    if (is.nan(SNARE_interactions_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],component])) { ##
      SNARE_interactions_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],component] <- sum(rna_file[rna_file$gene_name == component,"unstranded"])
    }
    else {
      cur <- SNARE_interactions_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],component]
      SNARE_interactions_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],component] <- round(cur + sum(rna_file[rna_file$gene_name == component,"unstranded"])/2)
    }
  }
}
#remove non-data gene
SNARE_interactions_frame <- SNARE_interactions_frame[,!is.nan(colSums(SNARE_interactions_frame))]
#feature scaling
process <- preProcess(SNARE_interactions_frame, method = c("range"))
scaled_SNARE_interactions_frame <- predict(process, SNARE_interactions_frame)

#####
Ubiquitin_frame <- data.frame(matrix(nrow = length(unique(RNA_sample_sheet[,"Case.ID"])), 
                                          ncol = length(Ubiquitin)))
rownames(Ubiquitin_frame) <- unique(RNA_sample_sheet[,"Case.ID"])
colnames(Ubiquitin_frame) <- Ubiquitin
Ubiquitin_frame[] <- NaN
pb <- progress_bar$new(total = length(All_RNA_dir)) #checking for loop progress
for (dir_name in All_RNA_dir) {
  pb$tick()
  rna_file_name <- list.files(dir_name, pattern = "*.tsv")
  if (!(rna_file_name %in% RNA_sample_sheet$File.Name)) {
    next
  }
  rna_file <- read.delim(paste0(dir_name,"/",rna_file_name),
                         sep = '\t',header = T,fill = T, skip = 1)
  for (component in Ubiquitin) {
    if (!(component %in% rna_file$gene_name)) {
      next
    }
    if (is.nan(Ubiquitin_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],component])) { ##
      Ubiquitin_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],component] <- sum(rna_file[rna_file$gene_name == component,"unstranded"])
    }
    else {
      cur <- Ubiquitin_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],component]
      Ubiquitin_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],component] <- round(cur + sum(rna_file[rna_file$gene_name == component,"unstranded"])/2)
    }
  }
}
#remove non-data gene
Ubiquitin_frame <- Ubiquitin_frame[,!is.nan(colSums(Ubiquitin_frame))]
#feature scaling
process <- preProcess(Ubiquitin_frame, method = c("range"))
scaled_Ubiquitin_frame <- predict(process, Ubiquitin_frame)

#####
Sulfur_relay_system_frame <- data.frame(matrix(nrow = length(unique(RNA_sample_sheet[,"Case.ID"])), 
                                          ncol = length(Sulfur_relay_system)))
rownames(Sulfur_relay_system_frame) <- unique(RNA_sample_sheet[,"Case.ID"])
colnames(Sulfur_relay_system_frame) <- Sulfur_relay_system
Sulfur_relay_system_frame[] <- NaN
pb <- progress_bar$new(total = length(All_RNA_dir)) #checking for loop progress
for (dir_name in All_RNA_dir) {
  pb$tick()
  rna_file_name <- list.files(dir_name, pattern = "*.tsv")
  if (!(rna_file_name %in% RNA_sample_sheet$File.Name)) {
    next
  }
  rna_file <- read.delim(paste0(dir_name,"/",rna_file_name),
                         sep = '\t',header = T,fill = T, skip = 1)
  for (component in Sulfur_relay_system) {
    if (!(component %in% rna_file$gene_name)) {
      next
    }
    if (is.nan(Sulfur_relay_system_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],component])) { ##
      Sulfur_relay_system_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],component] <- sum(rna_file[rna_file$gene_name == component,"unstranded"])
    }
    else {
      cur <- Sulfur_relay_system_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],component]
      Sulfur_relay_system_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],component] <- round(cur + sum(rna_file[rna_file$gene_name == component,"unstranded"])/2)
    }
  }
}
#remove non-data gene
Sulfur_relay_system_frame <- Sulfur_relay_system_frame[,!is.nan(colSums(Sulfur_relay_system_frame))]
#feature scaling
process <- preProcess(Sulfur_relay_system_frame, method = c("range"))
scaled_Sulfur_relay_system_frame <- predict(process, Sulfur_relay_system_frame)

#####
Proteasome_frame <- data.frame(matrix(nrow = length(unique(RNA_sample_sheet[,"Case.ID"])), 
                                          ncol = length(Proteasome)))
rownames(Proteasome_frame) <- unique(RNA_sample_sheet[,"Case.ID"])
colnames(Proteasome_frame) <- Proteasome
Proteasome_frame[] <- NaN
pb <- progress_bar$new(total = length(All_RNA_dir)) #checking for loop progress
for (dir_name in All_RNA_dir) {
  pb$tick()
  rna_file_name <- list.files(dir_name, pattern = "*.tsv")
  if (!(rna_file_name %in% RNA_sample_sheet$File.Name)) {
    next
  }
  rna_file <- read.delim(paste0(dir_name,"/",rna_file_name),
                         sep = '\t',header = T,fill = T, skip = 1)
  for (component in Proteasome) {
    if (!(component %in% rna_file$gene_name)) {
      next
    }
    if (is.nan(Proteasome_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],component])) { ##
      Proteasome_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],component] <- sum(rna_file[rna_file$gene_name == component,"unstranded"])
    }
    else {
      cur <- Proteasome_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],component]
      Proteasome_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],component] <- round(cur + sum(rna_file[rna_file$gene_name == component,"unstranded"])/2)
    }
  }
}
#remove non-data gene
Proteasome_frame <- Proteasome_frame[,!is.nan(colSums(Proteasome_frame))]
#feature scaling
process <- preProcess(Proteasome_frame, method = c("range"))
scaled_Proteasome_frame <- predict(process, Proteasome_frame)

#####
RNA_degradation_frame <- data.frame(matrix(nrow = length(unique(RNA_sample_sheet[,"Case.ID"])), 
                                          ncol = length(RNA_degradation)))
rownames(RNA_degradation_frame) <- unique(RNA_sample_sheet[,"Case.ID"])
colnames(RNA_degradation_frame) <- RNA_degradation
RNA_degradation_frame[] <- NaN
pb <- progress_bar$new(total = length(All_RNA_dir)) #checking for loop progress
for (dir_name in All_RNA_dir) {
  pb$tick()
  rna_file_name <- list.files(dir_name, pattern = "*.tsv")
  if (!(rna_file_name %in% RNA_sample_sheet$File.Name)) {
    next
  }
  rna_file <- read.delim(paste0(dir_name,"/",rna_file_name),
                         sep = '\t',header = T,fill = T, skip = 1)
  for (component in RNA_degradation) {
    if (!(component %in% rna_file$gene_name)) {
      next
    }
    if (is.nan(RNA_degradation_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],component])) { ##
      RNA_degradation_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],component] <- sum(rna_file[rna_file$gene_name == component,"unstranded"])
    }
    else {
      cur <- RNA_degradation_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],component]
      RNA_degradation_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],component] <- round(cur + sum(rna_file[rna_file$gene_name == component,"unstranded"])/2)
    }
  }
}
#remove non-data gene
RNA_degradation_frame <- RNA_degradation_frame[,!is.nan(colSums(RNA_degradation_frame))]
#feature scaling
process <- preProcess(RNA_degradation_frame, method = c("range"))
scaled_RNA_degradation_frame <- predict(process, RNA_degradation_frame)

standard_RNA_polymerase <- standardize(RNA_polymerase_frame)
standard_Basal_transcription_factors <- standardize(Basal_transcription_factors_frame)
standard_Spliceosome <- standardize(Spliceosome_frame)
standard_Ribosome <- standardize(Ribosome_frame)
standard_Aminoacyl_tRNA_biosynthesis <- standardize(Aminoacyl_tRNA_biosynthesis_frame)
standard_Nucleocytoplasmic_transport <- standardize(Nucleocytoplasmic_transport_frame)
standard_mRNA_surveillance_pathway <- standardize(mRNA_surveillance_pathway_frame)
standard_Ribosome_biogenesis_in_eukaryotes <- standardize(Ribosome_biogenesis_in_eukaryotes_frame)
standard_Protein_export <- standardize(Protein_export_frame)
standard_Protein_processing_in_endoplasmic_reticulum <- standardize(Protein_processing_in_endoplasmic_reticulum_frame)
standard_SNARE_interactions <- standardize(SNARE_interactions_frame)
standard_Ubiquitin <- standardize(Ubiquitin_frame)
standard_Sulfur_relay_system <- standardize(Sulfur_relay_system_frame)
standard_Proteasome <- standardize(Proteasome_frame)
standard_RNA_degradation <- standardize(RNA_degradation_frame)

level4 <- read.table(paste0(data_dir,'/Protein/TCGA-UCEC-L4.csv'), sep = ',', header = T)
level4 <- level4 %>% 
  filter(str_detect(Sample_ID, "-01A")) %>%
  separate(Sample_ID,c('Case_ID','dump'), sep = '-01A', remove = FALSE, extra = 'merge')
level4 <- level4[,c("Case_ID","ARID1A")]
level4 <- rename(level4, "Row.names" = "Case_ID")
level4 <- level4[!(level4$Row.names=="TCGA-EY-A1GJ"),]

scaled_RNA_polymerase_frame$Row.names <- rownames(scaled_RNA_polymerase_frame)
scaled_RNA_polymerase_frame <- merge(scaled_RNA_polymerase_frame, level4, by = 'Row.names', all.x = T)
scaled_Basal_transcription_factors_frame$Row.names <- rownames(scaled_Basal_transcription_factors_frame)
scaled_Basal_transcription_factors_frame <- merge(scaled_Basal_transcription_factors_frame, level4, by = 'Row.names', all.x = T)
scaled_Spliceosome_frame$Row.names <- rownames(scaled_Spliceosome_frame)
scaled_Spliceosome_frame <- merge(scaled_Spliceosome_frame, level4, by = 'Row.names', all.x = T)
scaled_Ribosome_frame$Row.names <- rownames(scaled_Ribosome_frame)
scaled_Ribosome_frame <- merge(scaled_Ribosome_frame, level4, by = 'Row.names', all.x = T)
scaled_Aminoacyl_tRNA_biosynthesis_frame$Row.names <- rownames(scaled_Aminoacyl_tRNA_biosynthesis_frame)
scaled_Aminoacyl_tRNA_biosynthesis_frame <- merge(scaled_Aminoacyl_tRNA_biosynthesis_frame, level4, by = 'Row.names', all.x = T)
scaled_Nucleocytoplasmic_transport_frame$Row.names <- rownames(scaled_Nucleocytoplasmic_transport_frame)
scaled_Nucleocytoplasmic_transport_frame <- merge(scaled_Nucleocytoplasmic_transport_frame, level4, by = 'Row.names', all.x = T)
scaled_mRNA_surveillance_pathway_frame$Row.names <- rownames(scaled_mRNA_surveillance_pathway_frame)
scaled_mRNA_surveillance_pathway_frame <- merge(scaled_mRNA_surveillance_pathway_frame, level4, by = 'Row.names', all.x = T)
scaled_Ribosome_biogenesis_in_eukaryotes_frame$Row.names <- rownames(scaled_Ribosome_biogenesis_in_eukaryotes_frame)
scaled_Ribosome_biogenesis_in_eukaryotes_frame <- merge(scaled_Ribosome_biogenesis_in_eukaryotes_frame, level4, by = 'Row.names', all.x = T)
scaled_Protein_export_frame$Row.names <- rownames(scaled_Protein_export_frame)
scaled_Protein_export_frame <- merge(scaled_Protein_export_frame, level4, by = 'Row.names', all.x = T)
scaled_Protein_processing_in_endoplasmic_reticulum_frame$Row.names <- rownames(scaled_Protein_processing_in_endoplasmic_reticulum_frame)
scaled_Protein_processing_in_endoplasmic_reticulum_frame <- merge(scaled_Protein_processing_in_endoplasmic_reticulum_frame, level4, by = 'Row.names', all.x = T)
scaled_SNARE_interactions_frame$Row.names <- rownames(scaled_SNARE_interactions_frame)
scaled_SNARE_interactions_frame <- merge(scaled_SNARE_interactions_frame, level4, by = 'Row.names', all.x = T)
scaled_Ubiquitin_frame$Row.names <- rownames(scaled_Ubiquitin_frame)
scaled_Ubiquitin_frame <- merge(scaled_Ubiquitin_frame, level4, by = 'Row.names', all.x = T)
scaled_Sulfur_relay_system_frame$Row.names <- rownames(scaled_Sulfur_relay_system_frame)
scaled_Sulfur_relay_system_frame <- merge(scaled_Sulfur_relay_system_frame, level4, by = 'Row.names', all.x = T)
scaled_Proteasome_frame$Row.names <- rownames(scaled_Proteasome_frame)
scaled_Proteasome_frame <- merge(scaled_Proteasome_frame, level4, by = 'Row.names', all.x = T)
scaled_RNA_degradation_frame$Row.names <- rownames(scaled_RNA_degradation_frame)
scaled_RNA_degradation_frame <- merge(scaled_RNA_degradation_frame, level4, by = 'Row.names', all.x = T)
# 
write.csv(scaled_RNA_polymerase_frame, paste0(data_dir,'/KEGG_RNA_polymerase_frame.csv'), sep = ',')
write.csv(scaled_Basal_transcription_factors_frame, paste0(data_dir,'/KEGG_Basal_transcription_factors_frame.csv'), sep = ',')
write.csv(scaled_Spliceosome_frame, paste0(data_dir,'/KEGG_Spliceosome_frame.csv'), sep = ',')
write.csv(scaled_Ribosome_frame, paste0(data_dir,'/KEGG_Ribosome_frame.csv'), sep = ',')
write.csv(scaled_Aminoacyl_tRNA_biosynthesis_frame, paste0(data_dir,'/KEGG_Aminoacyl_tRNA_biosynthesis_frame.csv'), sep = ',')
write.csv(scaled_Nucleocytoplasmic_transport_frame, paste0(data_dir,'/KEGG_Nucleocytoplasmic_transport_frame.csv'), sep = ',')
write.csv(scaled_mRNA_surveillance_pathway_frame, paste0(data_dir,'/KEGG_mRNA_surveillance_pathway_frame.csv'), sep = ',')
write.csv(scaled_Ribosome_biogenesis_in_eukaryotes_frame, paste0(data_dir,'/KEGG_Ribosome_biogenesis_in_eukaryotes_frame.csv'), sep = ',')
write.csv(scaled_Protein_export_frame, paste0(data_dir,'/KEGG_Protein_export_frame.csv'), sep = ',')
write.csv(scaled_Protein_processing_in_endoplasmic_reticulum_frame, paste0(data_dir,'/KEGG_Protein_processing_in_endoplasmic_reticulum_frame.csv'), sep = ',')
write.csv(scaled_SNARE_interactions_frame, paste0(data_dir,'/KEGG_SNARE_interactions_frame.csv'), sep = ',')
write.csv(scaled_Ubiquitin_frame,paste0(data_dir, '/KEGG_Ubiquitin_frame.csv'), sep = ',')
write.csv(scaled_Sulfur_relay_system_frame, paste0(data_dir,'/KEGG_Sulfur_relay_system_frame.csv'), sep = ',')
write.csv(scaled_Proteasome_frame, paste0(data_dir,'/KEGG_Proteasome_frame.csv'), sep = ',')
write.csv(scaled_RNA_degradation_frame, paste0(data_dir,'/KEGG_RNA_degradation_frame.csv'), sep = ',')


library("progress") 
library("tidyverse") 
library("dplyr")
library("TCGAbiolinks") 
library("stringr") 
library("TCGAutils") 
library('tibble') 
library('purrr')
library('caret') 
library(ggpubr)
library(cowplot)

data_dir <- "dataset"
RNA_sample_sheet_path <- paste0(data_dir, "/RNA_sample_sheet.tsv")
#A3 AO572
target_df <- readxl::read_xlsx(
  path = paste0(data_dir,'/Nature_2020_ARID1A_direct_targets.xlsx'),
  range = 'A3:AO572'
)
target_df <- target_df[3:569,c(1,9,19,29,39)]
target_df[,c(2,3,4,5)] <- target_df[,c(2,3,4,5)] > 0
target_df[,'Label'] <- rowSums(target_df[,c(2,3,4,5)])
target_df[,'Label'] <- target_df[,'Label'] > 2
target_df <- target_df[,c(1,6)]


#####################
## mRNA processing ##
#####################
RNA_sample_sheet <- read.delim(RNA_sample_sheet_path, sep = "\t")
Normal_sample_sheet <- RNA_sample_sheet %>% 
  filter(Sample.Type == 'Solid Tissue Normal')
RNA_sample_sheet <- RNA_sample_sheet %>% 
  filter(Sample.Type == "Primary Tumor") %>% 
  filter(str_detect(Sample.ID, "-01A"))

Patient_data_frame <- data.frame(matrix(nrow = length(unique(RNA_sample_sheet[,"Case.ID"])), ncol = length(target_df$`Gene Symbol`)))      ##filtering unique case ID
Normal_data_frame <- data.frame(matrix(nrow = length(unique(Normal_sample_sheet[,"Case.ID"])), ncol = length(target_df$`Gene Symbol`)))     ##filtering unique case ID

rownames(Patient_data_frame) <- unique(RNA_sample_sheet[,"Case.ID"])
colnames(Patient_data_frame) <- target_df$`Gene Symbol`
rownames(Normal_data_frame) <- unique(Normal_sample_sheet[,"Case.ID"])
colnames(Normal_data_frame) <- target_df$`Gene Symbol`

RNA_dir <- paste0(data_dir,"/RNA_seq")
All_RNA_dir <- list.dirs(RNA_dir, recursive = FALSE)

pb <- progress_bar$new(total = length(All_RNA_dir)) #checking for loop progress
for (dir_name in All_RNA_dir) {
  pb$tick()
  rna_file_name <- list.files(dir_name, pattern = "*.tsv")
  if (!(rna_file_name %in% RNA_sample_sheet$File.Name)) {
    next
  }
  rna_file <- read.delim(paste0(dir_name,"/",rna_file_name),
                         sep = '\t',header = T,fill = T, skip = 1)
  for (target in target_df$`Gene Symbol`) {
    if (!(target %in% rna_file$gene_name)) {
      next
    }
    target_count <- rna_file[rna_file$gene_name == target,"tpm_unstranded"]
    if (is.na(Patient_data_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],target])) { ##
      Patient_data_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],target] <- target_count
    }
    else {
      cur <- Patient_data_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],target]
      Patient_data_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],target] <- round((cur + target_count)/2)
    }
  }
}


pb <- progress_bar$new(total = length(All_RNA_dir)) #checking for loop progress
for (dir_name in All_RNA_dir) {
  pb$tick()
  rna_file_name <- list.files(dir_name, pattern = "*.tsv")
  if (!(rna_file_name %in% Normal_sample_sheet$File.Name)) {
    next
  }
  rna_file <- read.delim(paste0(dir_name,"/",rna_file_name),
                         sep = '\t',header = T,fill = T, skip = 1)
  for (target in target_df$`Gene Symbol`) {
    if (!(target %in% rna_file$gene_name)) {
      next
    }
    target_count <- rna_file[rna_file$gene_name == target,"tpm_unstranded"]
    if (is.na(Normal_data_frame[Normal_sample_sheet[Normal_sample_sheet$File.Name == rna_file_name,"Case.ID"],target])) { ##
      Normal_data_frame[Normal_sample_sheet[Normal_sample_sheet$File.Name == rna_file_name,"Case.ID"],target] <- target_count
    }
    else {
      cur <- Normal_data_frame[Normal_sample_sheet[Normal_sample_sheet$File.Name == rna_file_name,"Case.ID"],target]
      Normal_data_frame[Normal_sample_sheet[Normal_sample_sheet$File.Name == rna_file_name,"Case.ID"],target] <- round((cur + target_count)/2)
    }
  }
}

colnames(Patient_data_frame)[is.na(colSums(Patient_data_frame))]
colnames(Normal_data_frame)[is.na(colSums(Normal_data_frame))]

### remove NA genes
Patient_data_frame <- Patient_data_frame[,!(names(Patient_data_frame)%in%(colnames(Patient_data_frame)[is.na(colSums(Patient_data_frame))]))]
Normal_data_frame <- Normal_data_frame[,!(names(Normal_data_frame)%in%(colnames(Normal_data_frame)[is.na(colSums(Normal_data_frame))]))]

###############
## Labeling ###
###############
#### ARID1A Function Active = True
#### ARID1A Function Inactive = False
Normal_average <- colMeans(Normal_data_frame)
Patient_label <- data.frame(matrix(nrow = nrow(Patient_data_frame),ncol = ncol(Patient_data_frame)))
rownames(Patient_label) <- unique(RNA_sample_sheet[,"Case.ID"])
colnames(Patient_label) <- colnames(Patient_data_frame)
Patient_data_frame[1,] > Normal_average

UP_genes <- target_df$`Gene Symbol`[target_df$Label == TRUE]
UP_genes <- UP_genes[UP_genes %in% colnames(Patient_data_frame)]
DOWN_genes <- target_df$`Gene Symbol`[target_df$Label == FALSE]
DOWN_genes <- DOWN_genes[DOWN_genes %in% colnames(Patient_data_frame)]

## UP regulated gene labeling 
Patient_label[,UP_genes] <- Patient_data_frame[,UP_genes] < Normal_average[UP_genes]

## DOWN regulated gene labeling 
Patient_label[,DOWN_genes] <- Patient_data_frame[,DOWN_genes] > Normal_average[DOWN_genes]


Active_patients <- colSums(Patient_label)
hist(Active_patients, prob=TRUE, col="grey",breaks = 20)
lines(density(Active_patients), col="blue", lwd=2)

write.csv(Patient_label,file = paste0(data_dir,"/UCEC_554_targets_Functional_data.csv"))
write.csv(target_df,file = paste0(data_dir,"/554_targets_summary.csv"))

paste0(data_dir,"/RNA_seq")

# ARID1A extracting
Patient_data_frame$ARID1A <- NaN
pb <- progress_bar$new(total = length(All_RNA_dir)) #checking for loop progress
for (dir_name in All_RNA_dir) {
  pb$tick()
  rna_file_name <- list.files(dir_name, pattern = "*.tsv")
  if (!(rna_file_name %in% RNA_sample_sheet$File.Name)) {
    next
  }
  rna_file <- read.delim(paste0(dir_name,"/",rna_file_name),
                         sep = '\t',header = T,fill = T, skip = 1)
  ARID1A_count <- rna_file[rna_file$gene_name == "ARID1A","tpm_unstranded"]
  if (is.nan(Patient_data_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],"ARID1A"])) { ##
    Patient_data_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],"ARID1A"] <- ARID1A_count
  }
  else {
    cur <- Patient_data_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],"ARID1A"]
    Patient_data_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],"ARID1A"] <- round((cur + ARID1A_count)/2)
  }
}
Patient_data_frame['ARID1A']
write.csv(Patient_data_frame,file = paste0(data_dir,"/UCEC_mRNA_TPM_matrix.csv"))
#

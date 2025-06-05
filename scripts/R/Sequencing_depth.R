library("progress")
library("tidyverse") 
library("dplyr")
library("TCGAbiolinks") 
library("stringr") 
library("TCGAutils") 
library('tibble') 
library('caret') 
library(ggplot2)
library(ggpubr)
library(patchwork)


data_dir <- "dataset"
RNA_sample_sheet_path <- paste0(data_dir, "/RNA_sample_sheet.tsv")
RNA_sample_sheet <- read.delim(RNA_sample_sheet_path, sep = "\t")
RNA_sample_sheet <- RNA_sample_sheet %>% 
  filter(Sample.Type == "Primary Tumor") %>% 
  filter(str_detect(Sample.ID, "-01A"))
Patient_data_frame <- data.frame(row.names = unique(RNA_sample_sheet[,"Case.ID"]))      ##filtering unique case ID


#####################
## mRNA processing ##
#####################
Patient_data_frame$RNA_count <- NaN
Patient_data_frame$RNA_TPM <- NaN
Patient_data_frame$Seq_depth <- NaN
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
  ARID1A_count <- rna_file[rna_file$gene_name == "ARID1A","unstranded"]
  ARID1A_tpm <- rna_file[rna_file$gene_name == "ARID1A","tpm_unstranded"]
  seq_depth <- sum(rna_file$unstranded, na.rm = TRUE)
  if (is.nan(Patient_data_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],"RNA_count"])) { ##
    Patient_data_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],"RNA_count"] <- ARID1A_count
    Patient_data_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],"RNA_TPM"] <- ARID1A_tpm
    Patient_data_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],"Seq_depth"] <- seq_depth
  }
  else {
    cur <- Patient_data_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],"RNA_count"]
    cur_tpm <- Patient_data_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],"RNA_TPM"]
    cur_depth <- Patient_data_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],"Seq_depth"]
    Patient_data_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],"RNA_count"] <- round((cur + ARID1A_count)/2)
    Patient_data_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],"RNA_TPM"] <- round((cur + ARID1A_tpm)/2)
    Patient_data_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],"Seq_depth"] <- round((cur + seq_depth)/2)
  }
}

DEG_labels <- read.csv("dataset/UCEC_labeling_mRNA_protein_mutation_functional_state.csv")

cor.test(Patient_data_frame$RNA_count, Patient_data_frame$Seq_depth, method = "spearman")
cor.test(Patient_data_frame$RNA_TPM, Patient_data_frame$Seq_depth, method = "spearman")

## normalization
process <- preProcess(as.data.frame(Patient_data_frame$RNA_count), method = c("range"))
scaled_RNA <- predict(process, as.data.frame(Patient_data_frame$RNA_count))
Patient_data_frame$scaled_RNA_count <- scaled_RNA$`Patient_data_frame$RNA_count`

cor.test(Patient_data_frame$scaled_RNA_count, Patient_data_frame$Seq_depth, method = "spearman")


process <- preProcess(as.data.frame(Patient_data_frame$RNA_TPM), method = c("range"))
scaled_TPM <- predict(process, as.data.frame(Patient_data_frame$RNA_TPM))
Patient_data_frame$scaled_RNA_TPM <- scaled_TPM$`Patient_data_frame$RNA_TPM`

cor.test(Patient_data_frame$scaled_RNA_TPM, Patient_data_frame$scaled_RNA_count, method = "spearman")

Patient_data_frame$sample_id = rownames(Patient_data_frame)
Patient_data_frame <- Patient_data_frame %>% 
  left_join(DEG_labels, by= c('sample_id' = 'X'))

######
df_mrna <- Patient_data_frame %>% filter(mRNA_based_label %in% c("High", "Low"))
p1 <- ggplot(df_mrna, aes(x = mRNA_based_label, y = Seq_depth, fill = mRNA_based_label)) +
  geom_boxplot() +
  stat_compare_means(
    method = "kruskal.test",
    label = "p.signif",
    label.y = max(df_mrna$Seq_depth, na.rm = TRUE) * 0.95,
    label.x = 1.5,
    size = 7
  ) +
  labs(title = "mRNA-based Label", y = "Seq Depth", x = NULL) +
  theme_minimal(base_size = 14) +  
  theme(
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 15, face = "bold")
  )

# 2. Mutation Label plot
p2 <- ggplot(Patient_data_frame, aes(x = Mutation_label, y = Seq_depth, fill = Mutation_label)) +
  geom_boxplot() +
  stat_compare_means(
    method = "kruskal.test",
    label = "p.signif",
    label.y = max(Patient_data_frame$Seq_depth, na.rm = TRUE) * 0.95,
    label.x = 1.5,
    size = 7
  ) +
  labs(title = "Mutation Label", y = "Seq Depth", x = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 15, face = "bold")
  )

# 3. functional_state_based Label plot
df_functional_state_based <- Patient_data_frame %>% filter(functional_state_based_label %in% c("High", "Low"))
p3 <- ggplot(df_functional_state_based, aes(x = functional_state_based_label, y = Seq_depth, fill = functional_state_based_label)) +
  geom_boxplot() +
  stat_compare_means(
    method = "kruskal.test",
    label = "p.signif",
    label.y = max(df_functional_state_based$Seq_depth, na.rm = TRUE) * 0.95,
    label.x = 1.5,
    size = 7
  ) +
  labs(title = "Functional state based Label", y = "Seq Depth", x = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 15, face = "bold")
  )

# 4. Combine into one figure (horizontal)
combined_plot <- p1 + p2 + p3 + plot_layout(ncol = 3)

 


#####
## Compare sepuencing depth with predicted protein 
protein_df <- read.csv('dataset/UCEC_Observed_Predicted_ARID1A.csv')
Patient_data_frame <- Patient_data_frame %>% 
  left_join(protein_df, by= c('sample_id' = 'X'))
sum(is.na(Patient_data_frame$Observed_ARID1A))
sum(!is.na(Patient_data_frame$Observed_ARID1A))

Predicted_patients <- Patient_data_frame %>% filter(is.na(Observed_ARID1A))

cor.test(Predicted_patients$scaled_RNA_count, Predicted_patients$Predicted_ARID1A, method = "spearman")
cor.test(Predicted_patients$scaled_RNA_TPM, Predicted_patients$Predicted_ARID1A, method = "spearman")
cor.test(Predicted_patients$Seq_depth, Predicted_patients$Predicted_ARID1A, method = "spearman")
cor.test(Patient_data_frame$Seq_depth, Patient_data_frame$Predicted_ARID1A, method = "spearman")

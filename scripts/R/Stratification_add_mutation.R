library(DESeq2)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(cowplot)
library("tidyverse") 
library("dplyr")
library("TCGAbiolinks") 
library("stringr") 
library("BiocParallel")
library(biomaRt) 
library(org.Hs.eg.db)
library(AnnotationDbi)
register(SnowParam(6))
library(EnhancedVolcano)
library(factoextra)
library(cluster)

data_dir <- "dataset"

proj <- "TCGA-UCEC"
query <- GDCquery(
  project = proj,
  data.category = "Transcriptome Profiling", 
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
data <- GDCprepare(query)
DEG_labels <- read.csv(paste0(data_dir,"/UCEC_DEG_labels_mRNA_protein_functional_state_based_label.csv"))

# primary tumor sample filtering 
filtered_data <- data[,grep("-01A", colnames(data), value = TRUE)]
colnames(filtered_data) <- (sub("(-01A.*)", "", colnames(filtered_data)))
colnames(filtered_data) <- make.unique(colnames(filtered_data), sep = "__")
filtered_data <- filtered_data[,grep("__1", colnames(filtered_data), value = TRUE, invert = TRUE)]

# match the label order and Add ColData
rownames(DEG_labels) <- DEG_labels$X
DEG_labels <- DEG_labels[match(colnames(filtered_data), rownames(DEG_labels)),]
DEG_labels <- DEG_labels[,-1]
colData(filtered_data) <- cbind(colData(filtered_data), DEG_labels)

Mutation_label <- read.csv(paste0(data_dir,'/TCGA_UCEC_scaled.csv'))
# match the label order and Add ColData
rownames(Mutation_label) <- Mutation_label$X
Mutation_label <- Mutation_label[match(colnames(filtered_data), rownames(Mutation_label)),]
Mutation_label <- Mutation_label$Freq
Mutation_label <- if_else(Mutation_label > 0, "Mutated", "Not_mutated") 
colData(filtered_data) <- cbind(colData(filtered_data), Mutation_label)
#filtered_data$Mutation_label <- Mutation_label

############################
##### DATA EXPORT 
############################
Total_labels <- filtered_data@colData@listData[c(79:82)]
Total_labels <-  as.data.frame(Total_labels, row.names = colnames(filtered_data))
write.csv(Total_labels,paste0(data_dir,"/UCEC_labeling_mRNA_protein_mutation_functional_state.csv"))






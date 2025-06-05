library(ggplot2)
library(ggpubr)
library(ggsignif)
library(cowplot)
library("tidyverse") 
library("dplyr")
library(UpSetR)
library(dplyr)
library(TCGAbiolinks)
library(survminer)
library(survival)
library(SummarizedExperiment)
library(tidyverse)
library(gridExtra)
library(grid)


label_data <-  read.csv("dataset/UCEC_labeling_mRNA_protein_mutation_functional_state.csv")


MSI_df <- readxl::read_xlsx('dataset/TCGA_MSI_score.xlsx')
MSI_ucec <- MSI_df %>%
  dplyr::filter(`Cancer Type` == "TCGA-UCEC") %>%
  dplyr::select(`Case ID`, `MANTIS Score`)
MSI_ucec$Sample_ID <- MSI_ucec$`Case ID`  
label_data$Sample_ID <- label_data$X     
MSI_ucec <-  merge(MSI_ucec, label_data, all.x = TRUE, by='Sample_ID')
MSI_ucec <- MSI_ucec %>% 
  mutate(`MSI-H` = `MANTIS Score` > 0.4)

sum(MSI_ucec %>%
  filter(!is.na(Mutation_label)) %>% 
  dplyr::select(`MSI-H`))

MSI_ucec %>%
  group_by(mRNA_based_label) %>%
  summarise(MSI_H_ratio = mean(`MSI-H`, na.rm = TRUE))

MSI_ucec %>%
  group_by(functional_state_based_label) %>%
  summarise(MSI_H_ratio = mean(`MSI-H`, na.rm = TRUE))

MSI_ucec %>%
  group_by(Mutation_label) %>%
  summarise(MSI_H_ratio = mean(`MSI-H`, na.rm = TRUE))

MSI_ucec_filtered <- MSI_ucec %>%
  filter(
    (mRNA_based_label %in% c("High", "Low")) |
      (functional_state_based_label %in% c("High", "Low")) |
      (Mutation_label %in% c("Mutated", "Not_mutated"))
  )

MSI_summary <- bind_rows(
  MSI_ucec %>%
    filter(mRNA_based_label %in% c("High", "Low")) %>%
    group_by(Group = mRNA_based_label) %>%
    summarise(MSI_H_ratio = mean(`MSI-H`, na.rm = TRUE)) %>%
    mutate(Label_type = "mRNA_based_label"),
  
  MSI_ucec %>%
    filter(functional_state_based_label %in% c("High", "Low")) %>%
    group_by(Group = functional_state_based_label) %>%
    summarise(MSI_H_ratio = mean(`MSI-H`, na.rm = TRUE)) %>%
    mutate(Label_type = "functional_state_based_label"),
  
  MSI_ucec %>%
    filter(Mutation_label %in% c("Mutated", "Not_mutated")) %>%
    group_by(Group = Mutation_label) %>%
    summarise(MSI_H_ratio = mean(`MSI-H`, na.rm = TRUE)) %>%
    mutate(Label_type = "Mutation_label")
)


fisher.test(table(MSI_ucec$mRNA_based_label, MSI_ucec$`MSI-H`)[c("High", "Low"), ])
fisher.test(table(MSI_ucec$functional_state_based_label, MSI_ucec$`MSI-H`)[c("High", "Low"), ])
fisher.test(table(MSI_ucec$Mutation_label, MSI_ucec$`MSI-H`)[c("Mutated", "Not_mutated"), ])


signif_df <- data.frame(
  Label_type = c("mRNA_based_label", "Mutation_label", "functional_state_based_label"),
  xmin = c("High", "Mutated", "High"),
  xmax = c("Low", "Not_mutated", "Low"),
  annotations = c("NS", "***", "***"),
  y_position = c(0.41, 0.58, 0.61)
)

p <- ggplot(MSI_summary, aes(x = Group, y = MSI_H_ratio, fill = Group)) +
  geom_col(width = 1) +
  facet_wrap(~ Label_type, scales = "free_x") +
  theme_minimal(base_size = 14) +
  ylab("MSI-H Ratio") +
  ggtitle("MSI-H Ratio by Label Group") +
  scale_fill_brewer(palette = "Set2")

p + geom_signif(
  data = signif_df,
  aes(
    xmin = xmin,
    xmax = xmax,
    annotations = annotations,
    y_position = y_position
  ),
  manual = TRUE,
  inherit.aes = FALSE,
  textsize = 5,
  tip_length = 0.01
)

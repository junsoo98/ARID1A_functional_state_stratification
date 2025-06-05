library("progress")
library("tidyverse") 
library("dplyr")
library("TCGAbiolinks") 
library("stringr") 
library("TCGAutils") 
library('tibble') 
library('tidyr')
library('reshape2')
library(graphics)
library(ggplot2)
library(ggpubr)
library(openxlsx)
library(patchwork)


query <- GDCquery(
  project = "TCGA-UCEC", 
  data.category = "Simple Nucleotide Variation", 
  access = "open", 
  data.type = "Masked Somatic Mutation"
)
#GDCdownload(query)
SNV <- GDCprepare(query)
SNV_ARID1A <- SNV[SNV$Hugo_Symbol == "ARID1A",]
unique_samples <- unique(SNV$case_id)
sample_df <- data.frame(Sample = unique_samples, sample_id = unique_samples) 
for (UUID in unique_samples) {
  sample_df[sample_df$Sample == UUID,'sample_id'] <- UUIDtoBarcode(id_vector = UUID)$submitter_id
}

coding_size_mb <- 30

non_silent_classes <- c(
  "Missense_Mutation",
  "Nonsense_Mutation",
  "Frame_Shift_Del",
  "Frame_Shift_Ins",
  "Splice_Site",
  "In_Frame_Del",
  "In_Frame_Ins",
  "Nonstop_Mutation",
  "Translation_Start_Site"
)

tmb_df <- SNV %>%
  filter(Variant_Classification %in% non_silent_classes) %>%
  group_by(case_id) %>%
  summarise(
    Mutation_Count = n(),
    TMB = Mutation_Count / coding_size_mb
  ) %>%
  arrange(desc(TMB))

tmb_df <- tmb_df %>% 
  left_join(sample_df,by = c("case_id" = "Sample"))

DEG_labels <- read.csv("dataset/UCEC_labeling_mRNA_protein_mutation_functional_state.csv")

DEG_labels <- DEG_labels %>% 
  left_join(tmb_df, by= c('X' = 'sample_id')) %>% 
  filter(!is.na(TMB))


ggplot(DEG_labels %>% filter(mRNA_based_label %in% c("High", "Low")),
       aes(x = mRNA_based_label, y = TMB, fill = mRNA_based_label)) +
  stat_summary(fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.4) +
  stat_compare_means(
    comparisons = list(c("High", "Low")),
    label = "p.signif",
    method = "wilcox.test",
    label.y = 75  # 원하는 위치로 조정
  ) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(title = "Mean ± SE TMB by mRNA label", y = "TMB") +
  theme_minimal() + theme(legend.position = "none")

ggplot(DEG_labels%>% filter(Mutation_label %in% c("Not_mutated", "Mutated")),
       aes(x = Mutation_label, y = TMB, fill = Mutation_label)) +
  stat_summary(fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.4) +
  stat_compare_means(
    comparisons = list(c("Not_mutated", "Mutated")),
    label = "p.signif",
    method = "wilcox.test",
    label.y = 60  # 원하는 위치로 조정
  ) +  coord_cartesian(ylim = c(0, 100)) +
  labs(title = "Mean ± SE TMB by Mutation label", y = "TMB") +
  theme_minimal() + theme(legend.position = "none")

ggplot(DEG_labels %>% filter(functional_state_based_label %in% c("High", "Low")),
       aes(x = functional_state_based_label, y = TMB, fill = functional_state_based_label)) +
  stat_summary(fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  stat_compare_means(
    comparisons = list(c("High", "Low")),
    label = "p.signif",
    method = "wilcox.test",
    label.y = 75  # 원하는 위치로 조정
  ) +  coord_cartesian(ylim = c(0, 100)) +
  labs(title = "Mean ± SE TMB by functional_state_based label", y = "TMB") +
  theme_minimal() + theme(legend.position = "none")

####### combined 
# 1. mRNA-based
p1 <- ggplot(DEG_labels %>% filter(mRNA_based_label %in% c("High", "Low")),
             aes(x = mRNA_based_label, y = TMB, fill = mRNA_based_label)) +
  stat_summary(fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.4) +
  stat_compare_means(
    comparisons = list(c("High", "Low")),
    label = "p.signif",
    method = "wilcox.test",
    label.y = 60,
    size = 6
  ) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(title = "mRNA-based Label", y = "TMB", x = NULL) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

# 2. Mutation label
p2 <- ggplot(DEG_labels %>% filter(Mutation_label %in% c("Not_mutated", "Mutated")),
             aes(x = Mutation_label, y = TMB, fill = Mutation_label)) +
  stat_summary(fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.4) +
  stat_compare_means(
    comparisons = list(c("Not_mutated", "Mutated")),
    label = "p.signif",
    method = "wilcox.test",
    label.y = 60,
    size = 6
  ) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(title = "Mutation Label", y = "TMB", x = NULL) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

# 3. functional_state_based label
p3 <- ggplot(DEG_labels %>% filter(functional_state_based_label %in% c("High", "Low")),
             aes(x = functional_state_based_label, y = TMB, fill = functional_state_based_label)) +
  stat_summary(fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.4) +
  stat_compare_means(
    comparisons = list(c("High", "Low")),
    label = "p.signif",
    method = "wilcox.test",
    label.y = 60,
    size = 6
  ) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(title = "functional_state_based Label", y = "TMB", x = NULL) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

# 4. Combine horizontally
combined_plot <- p1 + p2 + p3 + plot_layout(ncol = 3)

combined_plot + plot_annotation(title = "TMB Comparison Across Labels")
######
summary_stats <- function(data, label_col) {
  data %>%
    group_by(!!sym(label_col)) %>%
    summarise(
      mean_TMB = mean(TMB, na.rm = TRUE),
      se_TMB = sd(TMB, na.rm = TRUE) / sqrt(sum(!is.na(TMB))),
      n = sum(!is.na(TMB))
    ) %>%
    mutate(label_type = label_col) %>%
    rename(label = !!sym(label_col))
}

# 각 레이블 컬럼별로 계산 후 결합
results <- bind_rows(
  summary_stats(DEG_labels, "mRNA_based_label"),
  summary_stats(DEG_labels, "functional_state_based_label"),
  summary_stats(DEG_labels, "Mutation_label")
)

# 보기 좋게 정렬
results <- results %>%
  select(label_type, label, mean_TMB, se_TMB, n)


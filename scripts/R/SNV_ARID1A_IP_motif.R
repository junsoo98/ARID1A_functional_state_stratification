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
library(SummarizedExperiment)



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

DEG_labels <- read.csv("dataset/UCEC_labeling_mRNA_protein_mutation_functional_state.csv")
rownames(DEG_labels) <- DEG_labels$X
protein_df <- read.csv('dataset/UCEC_Observed_Predicted_ARID1A.csv')
rownames(protein_df) <- protein_df$X
feature_df <- read.csv('dataset/TCGA_UCEC_scaled.csv')
rownames(feature_df) <- feature_df$X

observed_group <-  protein_df %>% filter(!is.na(Observed_ARID1A)) %>% select(X)

protein_SNV_df <- protein_df %>% left_join(feature_df[c("X","Freq","FS_SG","IMPACT","Severe_Position")], 
                         by = 'X')
cor.test(protein_SNV_df$Freq, protein_SNV_df$Predicted_ARID1A, method = "spearman")
cor.test(protein_SNV_df$FS_SG, protein_SNV_df$Predicted_ARID1A, method = "spearman")
cor.test(protein_SNV_df$IMPACT, protein_SNV_df$Predicted_ARID1A, method = "spearman")
cor.test(protein_SNV_df$Severe_Position, protein_SNV_df$Predicted_ARID1A, method = "spearman")

protein_SNV_df <- na.omit(protein_SNV_df)
cor.test(protein_SNV_df$Freq, protein_SNV_df$Observed_ARID1A, method = "spearman")
cor.test(protein_SNV_df$FS_SG, protein_SNV_df$Observed_ARID1A, method = "spearman")
cor.test(protein_SNV_df$IMPACT, protein_SNV_df$Observed_ARID1A, method = "spearman")
cor.test(protein_SNV_df$Severe_Position, protein_SNV_df$Observed_ARID1A, method = "spearman")

label_SNV_df <- DEG_labels %>% left_join(feature_df[c("X","Freq","FS_SG","IMPACT","Severe_Position")], 
                         by = 'X')
wilcox.test(Freq ~ functional_state_based_label, data = label_SNV_df %>% filter(functional_state_based_label != ""))
wilcox.test(FS_SG ~ functional_state_based_label, data = label_SNV_df %>% filter(functional_state_based_label != ""))
wilcox.test(IMPACT ~ functional_state_based_label, data = label_SNV_df %>% filter(functional_state_based_label != ""))
wilcox.test(Severe_Position ~ functional_state_based_label, data = label_SNV_df %>% filter(functional_state_based_label != ""))

wilcox.test(Freq ~ protein_based_label , data = label_SNV_df %>% filter(protein_based_label  %in% c('High','Low')))
wilcox.test(FS_SG ~ protein_based_label , data = label_SNV_df %>% filter(protein_based_label  %in% c('High','Low')))
wilcox.test(IMPACT ~ protein_based_label , data = label_SNV_df %>% filter(protein_based_label  %in% c('High','Low')))
wilcox.test(Severe_Position ~ protein_based_label , data = label_SNV_df %>% filter(protein_based_label  %in% c('High','Low')))

############################################
########### ARID1A IP binding motif 
############################################
library(GenomicFeatures)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(rtracklayer)

#ARID1A transcript CDS
txdb <- makeTxDbFromGFF("gencode.v43.annotation.gtf.gz", format = "gtf")
cds_list <- cdsBy(txdb, by = "tx", use.names = TRUE)
ARID1A_cds <- cds_list[['ENST00000324856.13']]

#HPA005456 epitope
ARID1A_protein <- readAAStringSet('dataset/O14497.fasta')
full_protein_seq <- ARID1A_protein[[1]]
epitope <- "PGLGNVAMGPRQHYPYGGPYDRVRTEPGIGPEGNMSTGAPQPNLMPSNPDSGMYSPSRYPPQQQQQQQQRHDSYGNQFSTQGTPSGSPFPSQQTTMYQQQQQNYK"

## epitope amino acid location
start_aa <- regexpr(epitope, full_protein_seq)[1]
end_aa <- start_aa + nchar(epitope)

# epitope nucleotide location 
cds_start_nt <- (start_aa - 1) * 3 + 1
cds_end_nt <- (end_aa) * 3

# epitope nt CDS exon mapping 
cds_ranges <- ARID1A_cds[order(start(ARID1A_cds))]
cds_width <- width(cds_ranges)
cds_cumsum <- cumsum(cds_width)
genomic_coords <- GRanges()
current_start <- cds_start_nt
current_end <- cds_end_nt
for (i in seq_along(cds_ranges)) {
  exon_start_nt <- ifelse(i == 1, 1, cds_cumsum[i - 1] + 1)
  exon_end_nt <- cds_cumsum[i]
  
  if (current_start <= exon_end_nt && current_end >= exon_start_nt) {
    overlap_start <- max(current_start, exon_start_nt)
    overlap_end <- min(current_end, exon_end_nt)
    offset_start <- overlap_start - exon_start_nt
    offset_end <- overlap_end - exon_start_nt
    
    exon_range <- cds_ranges[i]
    if (as.character(strand(exon_range)) == "+") {
      exon_genomic_start <- start(exon_range) + offset_start
      exon_genomic_end <- start(exon_range) + offset_end
    } else {
      exon_genomic_end <- end(exon_range) - offset_start
      exon_genomic_start <- end(exon_range) - offset_end
    }
    
    genomic_coords <- c(genomic_coords, GRanges(seqnames = seqnames(exon_range),
                                                ranges = IRanges(start = exon_genomic_start,
                                                                 end = exon_genomic_end),
                                                strand = strand(exon_range)))
  }
}
epitope_exons <- genomic_coords

#ARID1A SNV processing_only non silent 
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

SNV_ARID1A_nonsilient <- SNV_ARID1A %>%
  filter(Variant_Classification %in% non_silent_classes)
SNV_ranges <- GRanges(
  seqnames = SNV_ARID1A_nonsilient$Chromosome,
  ranges = IRanges(start = SNV_ARID1A_nonsilient$Start_Position,
                   end = SNV_ARID1A_nonsilient$End_Position),
  strand = SNV_ARID1A_nonsilient$Strand
)

overlap_hits <- findOverlaps(SNV_ranges, epitope_exons)
snvs_in_epitope <- SNV_ranges[queryHits(overlap_hits)]
SNV_ARID1A_in_epitope <- SNV_ARID1A_nonsilient[queryHits(overlap_hits), ]
### epitope around 50bp SNV
epitope_window <- resize(epitope_exons,
                         width = width(epitope_exons) + 100,
                         fix = "center")  # ±50bp

nearby_hits <- findOverlaps(SNV_ranges, epitope_window)
snvs_near_epitope <- SNV_ranges[queryHits(nearby_hits)]
SNV_ARID1A_near_epitope <- SNV_ARID1A_nonsilient[queryHits(nearby_hits), ]
## Distance from epitope center 
epitope_center <- mean(c(start(ranges(epitope_exons)), end(ranges(epitope_exons))))
snv_positions <- SNV_ARID1A_nonsilient$Start_Position
SNV_ARID1A_nonsilient$relative_to_epitope <- snv_positions - epitope_center

epitope_df <- as.data.frame(epitope_exons)
epitope_df$start_rel <- epitope_df$start - epitope_center
epitope_df$end_rel <- epitope_df$end - epitope_center


### Visualization 
library(ggplot2)
library(dplyr)

df_plot <- SNV_ARID1A_nonsilient %>%
  filter(abs(relative_to_epitope) <= 1000)  # ±1kb 내
#### Distance from epitope center around 1000bp
ggplot(df_plot, aes(x = relative_to_epitope)) +
  # 빨간 dashed 선: epitope 중심
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  
  # 회색 박스: epitope exon의 상대 위치 (epitope_df에서 직접 사용)
  geom_rect(data = epitope_df,
            inherit.aes = FALSE,  # 중요: df_plot에서 aes 상속하지 않음
            aes(xmin = start_rel, xmax = end_rel, ymin = 0, ymax = 0.4),
            fill = "gray70", alpha = 0.3) +
  
  # SNV 위치 점 (상대 거리)
  geom_point(aes(y = 0.2), position = position_jitter(height = 0.1), alpha = 0.7, color = "steelblue") +
  
  theme_minimal() +
  labs(x = "Distance from epitope center (bp)",
       y = NULL,
       title = "SNVs and epitope exon positions (±1000bp)") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

epitope_df <- as.data.frame(epitope_exons)
epitope_df$type <- "Epitope Exon"

## Epitope exon and SNV distribution in ARID1A gene 
snv_df <- SNV_ARID1A_nonsilient %>%
  mutate(pos = Start_Position,
         type = "SNV") %>%
  dplyr::select(pos, type)

viz_df <- rbind(
  data.frame(pos = (epitope_df$start + epitope_df$end) / 2, type = epitope_df$type),
  snv_df
)

ggplot(viz_df, aes(x = pos, y = type)) +
  geom_point(aes(color = type), size = 2, alpha = 0.8) +
  theme_minimal() +
  labs(x = "Genomic coordinate (chr1)", y = "",
       title = "SNVs and epitope exon locations") +
  scale_color_manual(values = c("Epitope Exon" = "red", "SNV" = "blue"))


# 1. SNV 위치 데이터 준비
snv_positions <- SNV_ARID1A_nonsilient$Start_Position

# 2. epitope exon을 표시하기 위한 데이터프레임
epitope_df <- as.data.frame(epitope_exons)

# 3. 히스토그램 + epitope 영역 오버레이
ggplot() +
  # 히스토그램 (빈도 기준)
  geom_histogram(aes(x = snv_positions, y = ..density..),
                 binwidth = 300,
                 fill = "skyblue", color = "black", alpha = 0.7) +
  
  # 밀도 곡선 (부드러운 분포)
  geom_density(aes(x = snv_positions), color = "darkblue", size = 1) +
  
  # epitope exon 위치 표시 (회색 박스)
  geom_rect(data = epitope_df,
            aes(xmin = start, xmax = end, ymin = 0, ymax = Inf),
            fill = "red", alpha = 0.3) +
  
  theme_minimal() +
  labs(title = "SNV density around ARID1A region; Red: epitope exon",
       x = "Genomic coordinate (chr1)",
       y = "Density of SNVs") +
  scale_x_continuous(labels = scales::comma)

### Protein expression who have epitome SNV  (vs whole sample )
SNV_ARID1A_near_epitope
snvs_near_epitope
sample_ids_50bp <- SNV_ARID1A[queryHits(nearby_hits), "Tumor_Sample_Barcode", drop = TRUE]
sample_ids_50bp <- str_remove(sample_ids_50bp, "-01A.*")
protein_df <- protein_df %>% mutate(group = ifelse(X %in% sample_ids_50bp, "With_epitope_SNV", "No_epitope_SNV" ))
table(protein_df %>% filter(Observed_ARID1A != 'NA') %>% select(group))

protein_df_na <- protein_df %>% filter(Observed_ARID1A != 'NA')

wilcox.test(Observed_ARID1A ~ group, data = protein_df_na)

obs_stat <- wilcox.test(Observed_ARID1A ~ group, data = protein_df_na)$statistic

perm_stats <- replicate(10000, {
  shuffled <- sample(protein_df_na$group)
  wilcox.test(protein_df_na$Observed_ARID1A ~ shuffled)$statistic
})

p_perm <- mean(perm_stats >= obs_stat)

### Protein expression who have epitome SNV (vs non silent SNV sample)
sample_ids_non_silent <- UUIDtoBarcode(unique(SNV_ARID1A_nonsilient$case_id))$submitter_id

protein_df_na_non_silent <- protein_df %>%
  filter(Observed_ARID1A != 'NA') %>% 
  filter(X %in% sample_ids_non_silent)
  
table(protein_df_na_non_silent$group)

wilcox.test(Observed_ARID1A ~ group, data = protein_df_na_non_silent)

obs_stat <- wilcox.test(Observed_ARID1A ~ group, data = protein_df_na_non_silent)$statistic

perm_stats <- replicate(10000, {
  shuffled <- sample(protein_df_na_non_silent$group)
  wilcox.test(protein_df_na_non_silent$Observed_ARID1A ~ shuffled)$statistic
})

p_perm <- mean(perm_stats >= obs_stat)


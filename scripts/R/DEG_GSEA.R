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
library("stringr") 
library("TCGAutils") 
library('tibble') 
library('caret') 

proj <- "TCGA-UCEC"
query <- GDCquery(
  project = proj,
  data.category = "Transcriptome Profiling", 
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
#GDCdownload(query)
data <- GDCprepare(query)
DEG_labels <- read.csv("dataset/UCEC_labeling_mRNA_protein_mutation_functional_state.csv")

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

#############################
######  1. mRNA
#############################

#filter NA annotation 
mRNA_data <- filtered_data[,!is.na(filtered_data$mRNA_based_label)]
mRNA_data <- mRNA_data[,!(filtered_data$mRNA_based_label == "Medium")] # Filter only High and Low

# make DESeq object
mRNA_ddsSE <- DESeqDataSet(mRNA_data, design = ~ mRNA_based_label)

# pre filtering low count genes (at least 10 counts, included in 3 patients)
smallestGroupSize <- 3
keep <- rowSums(counts(mRNA_ddsSE) >= 10) >= smallestGroupSize
mRNA_ddsSE <- mRNA_ddsSE[keep,]

# Apply DE analysis
mRNA_ddsSE <- DESeq(mRNA_ddsSE,parallel=TRUE)

# Take the result 
mRNA_res <- results(mRNA_ddsSE,parallel=TRUE)

# re-order by p-adj
mRNA_res <- mRNA_res[order(mRNA_res$padj),]

mRNA_gene_symbols <- mapIds(org.Hs.eg.db, sub("[.][0-9]*","",rownames(mRNA_res)), keytype = "ENSEMBL", column = "SYMBOL")
mRNA_res$ENSEMBL <- sub("[.][0-9]*","",rownames(mRNA_res))
mRNA_res$SYMBOL <- mRNA_gene_symbols


#summary 
mRNA_sig <- subset(mRNA_res, padj<.01 & abs(log2FoldChange)>1)

#############################
######  2. Protein
#############################

#filter NA annotation 
protein_data <- filtered_data[,!(filtered_data$protein_based_label == "nan")]
protein_data <- protein_data[,!(protein_data$protein_based_label == "Medium")] # Filter only High and Low

# make DESeq object
protein_ddsSE <- DESeqDataSet(protein_data, design = ~ protein_based_label)

# pre filtering low count genes (at least 10 counts, included in 3 patients)
smallestGroupSize <- 3
keep <- rowSums(counts(protein_ddsSE) >= 10) >= smallestGroupSize
protein_ddsSE <- protein_ddsSE[keep,]

# Apply DE analysis
protein_ddsSE <- DESeq(protein_ddsSE,parallel=TRUE)

# Take the result 
protein_res <- results(protein_ddsSE,parallel=TRUE)


# re-order by p-adj
protein_res <- protein_res[order(protein_res$padj),]

protein_gene_symbols <- mapIds(org.Hs.eg.db, sub("[.][0-9]*","",rownames(protein_res)), keytype = "ENSEMBL", column = "SYMBOL")
protein_res$ENSEMBL <- sub("[.][0-9]*","",rownames(protein_res))
protein_res$SYMBOL <- protein_gene_symbols


#summary 
protein_sig <- subset(protein_res, padj<.01 & abs(log2FoldChange)>1)


#############################
######  3. functional_state_based
#############################

#filter NA annotation 
functional_state_based_data <- filtered_data[,!is.na(filtered_data$functional_state_based_label)]
functional_state_based_data <- functional_state_based_data[,!(filtered_data$functional_state_based_label == "")] # Filter only High and Low

# make DESeq object
functional_state_based_ddsSE <- DESeqDataSet(functional_state_based_data, design = ~ functional_state_based_label)

# pre filtering low count genes (at least 10 counts, included in 3 patients)
smallestGroupSize <- 3
keep <- rowSums(counts(functional_state_based_ddsSE) >= 10) >= smallestGroupSize
functional_state_based_ddsSE <- functional_state_based_ddsSE[keep,]

# Apply DE analysis
functional_state_based_ddsSE <- DESeq(functional_state_based_ddsSE,parallel=TRUE)

# Take the result 
functional_state_based_res <- results(functional_state_based_ddsSE,parallel=TRUE)

# re-order by p-adj
functional_state_based_res <- functional_state_based_res[order(functional_state_based_res$padj),]


functional_state_based_gene_symbols <- mapIds(org.Hs.eg.db, sub("[.][0-9]*","",rownames(functional_state_based_res)), keytype = "ENSEMBL", column = "SYMBOL")
functional_state_based_res$ENSEMBL <- sub("[.][0-9]*","",rownames(functional_state_based_res))
functional_state_based_res$SYMBOL <- functional_state_based_gene_symbols

#summary 
functional_state_based_sig <- subset(functional_state_based_res, padj<.01 & abs(log2FoldChange)>1)


#############################
######  4. Mutation
#############################

#filter NA annotation 
Mutation_data <- filtered_data[,!is.na(filtered_data$Mutation_label)]

# make DESeq object
Mutation_ddsSE <- DESeqDataSet(Mutation_data, design = ~ Mutation_label)

# pre filtering low count genes (at least 10 counts, included in 3 patients)
smallestGroupSize <- 3
keep <- rowSums(counts(Mutation_ddsSE) >= 10) >= smallestGroupSize
Mutation_ddsSE <- Mutation_ddsSE[keep,]

# Apply DE analysis
Mutation_ddsSE <- DESeq(Mutation_ddsSE,parallel=TRUE)

# Take the result 
Mutation_res <- results(Mutation_ddsSE,parallel=TRUE)

# re-order by p-adj
Mutation_res <- Mutation_res[order(Mutation_res$padj),]


Mutation_gene_symbols <- mapIds(org.Hs.eg.db, sub("[.][0-9]*","",rownames(Mutation_res)), keytype = "ENSEMBL", column = "SYMBOL")
Mutation_res$ENSEMBL <- sub("[.][0-9]*","",rownames(Mutation_res))
Mutation_res$SYMBOL <- Mutation_gene_symbols

#summary 
Mutation_sig <- subset(Mutation_res, padj<.01 & abs(log2FoldChange)>1)

#############################
#### DEG visualization
#############################
DEG_list <- list(
  mRNA_label = list(rownames(mRNA_sig))[[1]],
  Mutation_label = list(rownames(Mutation_sig))[[1]],
  functional_state_label = list(rownames(functional_state_based_sig))[[1]]
)
display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL,...)
  grid.draw(venn_object)
}
display_venn(DEG_list,fill = c("#E69F00", "#56B4E9",'pink'),cat.cex = 1.4, cex = 2,
             fontfamily ="serif",main="Relationship between Significant DEGs", main.fontfamily="serif", main.cex=2.5)

#############
### GSEA & ORA
#############

########
## mRNA
########
gsea_res <- mRNA_res[order(-mRNA_res$stat),]
gene_list <- gsea_res$stat
names(gene_list) <- gsea_res$ENSEMBL
mRNA_gse <- gseGO(gene_list,
                  ont ="BP", 
                  keyType = "ENSEMBL", 
                  minGSSize = 10, 
                  maxGSSize = 500, 
                  pvalueCutoff = 0.05, 
                  OrgDb = 'org.Hs.eg.db'
)

gseaplot(mRNA_gse, geneSetID = 150)
mRNA_gse_term <- pairwise_termsim(mRNA_gse)
emapplot(mRNA_gse_term, showCategory=20, color = 'enrichmentScore', layout = 'kk',cex_label_category = 1.25)
emapplot_cluster(mRNA_gse_term, showCategory=20, color = 'enrichmentScore', layout = 'kk',cex_label_group = 1)
dotplot(mRNA_gse, showCategory=20, title = 'mRNA GSEA dotplot', font.size = 15, label_format = 50)


mRNA_sig_gsea_res <- mRNA_sig[order(-mRNA_sig$stat),]
mRNA_sig_gene_list <- mRNA_sig_gsea_res$stat
names(mRNA_sig_gene_list) <- mRNA_sig_gsea_res$ENSEMBL
mRNA_ggo <- groupGO(gene  = names(mRNA_sig_gene_list),
                    OrgDb    = org.Hs.eg.db,
                    keyType = 'ENSEMBL',
                    ont      = "BP",
                    level    = 3,
                    readable = TRUE)


mRNA_ego <- enrichGO(gene = names(mRNA_sig_gene_list),
                     universe = names(gene_list),
                     keyType = "ENSEMBL",
                     OrgDb = org.Hs.eg.db,
                     ont = "BP",
                     pAdjustMethod = "BH",
                     qvalueCutoff = 0.05,
                     readable = TRUE)

dotplot(mRNA_ego, showCategory=20, title = 'mRNA significant DEGs', font.size = 15, label_format = 50)
barplot(mRNA_ego, showCategory=20, title = 'mRNA significant DEGs', font.size = 15, label_format = 50)
mRNA_ego_term <- pairwise_termsim(mRNA_ego)
emapplot(mRNA_ego_term, showCategory=20,  layout = 'kk',cex_label_category = 1.5)
emapplot_cluster(mRNA_ego_term, showCategory=20, layout = 'kk', cex_label_group = 2.5)
mRNA_signif_res_lFC <- mRNA_sig_gsea_res$log2FoldChange
cnetplot(mRNA_ego,
         showCategory = 5,
         foldChange= mRNA_signif_res_lFC,
         vertex.label.font=6,
         circular = TRUE,
         colorEdge = TRUE
)

plotMA(mRNA_res, ylim=c(-5,5))


########
## functional_state_based
########
gsea_res <- functional_state_based_res[order(-functional_state_based_res$stat),]
gene_list <- gsea_res$stat
names(gene_list) <- gsea_res$ENSEMBL

functional_state_based_gse <- gseGO(gene_list,
                   ont ="BP", 
                   keyType = "ENSEMBL", 
                   minGSSize = 10, 
                   maxGSSize = 500, 
                   pvalueCutoff = 0.05, 
                   OrgDb = 'org.Hs.eg.db',
)

gseaplot(functional_state_based_gse, geneSetID = 150)
functional_state_based_gse_term <- pairwise_termsim(functional_state_based_gse)
emapplot(functional_state_based_gse_term, showCategory=20, color = 'enrichmentScore', layout = 'kk',cex_label_category = 1.25)
emapplot_cluster(functional_state_based_gse_term, showCategory=20, color = 'enrichmentScore', layout = 'kk',cex_label_group = 1)
dotplot(functional_state_based_gse, showCategory=20, title = 'functional_state_based GSEA dotplot', font.size = 15, label_format = 50)


functional_state_based_sig_gsea_res <- functional_state_based_sig[order(-functional_state_based_sig$stat),]
functional_state_based_sig_gene_list <- functional_state_based_sig_gsea_res$stat
names(functional_state_based_sig_gene_list) <- functional_state_based_sig_gsea_res$ENSEMBL

functional_state_based_ggo <- groupGO(gene     = names(functional_state_based_sig_gene_list),
                     OrgDb    = org.Hs.eg.db,
                     keyType = 'ENSEMBL',
                     ont      = "BP",
                     level    = 3,
                     readable = TRUE)

functional_state_based_ego <- enrichGO(gene = names(functional_state_based_sig_gene_list),
                      universe = names(gene_list),
                      keyType = "ENSEMBL",
                      OrgDb = org.Hs.eg.db,
                      ont = "BP",
                      pAdjustMethod = "BH",
                      qvalueCutoff = 0.05,
                      readable = TRUE)
dotplot(functional_state_based_ego, showCategory=20, title = 'functional state method significant DEGs', font.size = 15, label_format = 50)
barplot(functional_state_based_ego, showCategory=20, title = 'functional state method significant DEGs', font.size = 15, label_format = 50, )
functional_state_based_ego_term <- pairwise_termsim(functional_state_based_ego)
emapplot(functional_state_based_ego_term, showCategory=20, layout = 'kk',cex_label_category = 1.5)
emapplot_cluster(functional_state_based_ego_term, showCategory=20,layout = 'kk',cex_label_group = 2.5)
#emapplot(functional_state_based_ego, showCategory=10)
functional_state_based_signif_res_lFC <- functional_state_based_sig_gsea_res$log2FoldChange
cnetplot(functional_state_based_ego,
         showCategory = 5,
         foldChange = functional_state_based_signif_res_lFC,
         circular = TRUE,
         colorEdge = TRUE
)

plotMA(functional_state_based_res, ylim=c(-5,5))

###################
########
## Mutation
########
gsea_res <- Mutation_res[order(-Mutation_res$stat),]
gene_list <- gsea_res$stat
names(gene_list) <- gsea_res$ENSEMBL

Mutation_gse <- gseGO(gene_list,
                      ont ="BP", 
                      keyType = "ENSEMBL", 
                      minGSSize = 10, 
                      maxGSSize = 500, 
                      pvalueCutoff = 0.05, 
                      OrgDb = 'org.Hs.eg.db',
)

gseaplot(Mutation_gse, geneSetID = 150)
Mutation_gse_term <- pairwise_termsim(Mutation_gse)
emapplot(Mutation_gse_term, showCategory=20, color = 'enrichmentScore', layout = 'kk',cex_label_category = 1.25)
emapplot_cluster(Mutation_gse_term, showCategory=20, color = 'enrichmentScore', layout = 'kk',cex_label_group = 1)
dotplot(Mutation_gse, showCategory=20, title = 'Mutation GSEA dotplot', font.size = 15, label_format = 50)


Mutation_sig_gsea_res <- Mutation_sig[order(-Mutation_sig$stat),]
Mutation_sig_gene_list <- Mutation_sig_gsea_res$stat
names(Mutation_sig_gene_list) <- Mutation_sig_gsea_res$ENSEMBL

Mutation_ggo <- groupGO(gene     = names(Mutation_sig_gene_list),
                        OrgDb    = org.Hs.eg.db,
                        keyType = 'ENSEMBL',
                        ont      = "BP",
                        level    = 3,
                        readable = TRUE)

Mutation_ego <- enrichGO(gene = names(Mutation_sig_gene_list),
                         universe = names(gene_list),
                         keyType = "ENSEMBL",
                         OrgDb = org.Hs.eg.db,
                         ont = "BP",
                         pAdjustMethod = "BH",
                         qvalueCutoff = 0.05,
                         readable = TRUE)
dotplot(Mutation_ego, showCategory=20,title = 'Mutation label significant DEGs',font.size = 15, label_format = 50)
barplot(Mutation_ego, showCategory=20, title = 'Mutation label significant DEGs', font.size = 15, label_format = 50 )
Mutation_ego_term <- pairwise_termsim(Mutation_ego)
emapplot(Mutation_ego_term, showCategory=20, layout = 'kk',cex_label_category = 1.5)
emapplot_cluster(Mutation_ego_term, showCategory=20,layout = 'kk',cex_label_group = 2.5)
#emapplot(Mutation_ego, showCategory=10)
Mutation_signif_res_lFC <- Mutation_sig_gsea_res$log2FoldChange
cnetplot(Mutation_ego,
         showCategory = 5,
         foldChange = Mutation_signif_res_lFC,
         circular = TRUE,
         colorEdge = TRUE
)

plotMA(Mutation_res, ylim=c(-5,5))


plot_gene_box <- function(dds_object, gene_list, group_col = "protein_based_label", High = "High", Low = 'Low') {
  plots <- list()
  
  for (i in seq_along(gene_list)) {
    gene_name <- names(gene_list)[i]
    gene_id <- gene_list[[i]]
    
    df <- plotCounts(dds_object, gene = gene_id, intgroup = group_col, returnData = TRUE)
    
    p <- ggplot(df, aes_string(x = group_col, y = "count")) +
      geom_boxplot() +
      geom_signif(comparisons = list(c(High, Low)), map_signif_level = TRUE) +
      ggtitle(gene_name)
    
    plots[[i]] <- p
  }
  
  return(plot_grid(plotlist = plots, labels = "AUTO"))
}

gene_list <- list(
  ARID1A = "ENSG00000117713.20",
  TUBB = "ENSG00000196230.13",
  GATA3 = "ENSG00000107485.18",
  CD274 = "ENSG00000120217.14",
  ESR1 = "ENSG00000091831.24",
  FOXM1 = "ENSG00000111206.13",
  CXCL9 = "ENSG00000138755.6",
  PIK3CA = "ENSG00000121879.6",
  FOXA1 = "ENSG00000129514.8",
  KRAS = "ENSG00000133703.13",
  HLA_A = "ENSG00000206503.13",
  HLA_B = "ENSG00000234745.12",
  EZH2 = 'ENSG00000106462.11',
  PTEN = 'ENSG00000171862.11',
  TP53 = 'ENSG00000141510.18'
)

plot_gene_box(mRNA_ddsSE, gene_list, group_col = "mRNA_based_label")
plot_gene_box(Mutation_ddsSE, gene_list, group_col = "Mutation_label", High = "Not_mutated", Low = "Mutated")
plot_gene_box(functional_state_based_ddsSE, gene_list, group_col = "functional_state_based_label")
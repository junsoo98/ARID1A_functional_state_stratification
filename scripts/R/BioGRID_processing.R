library("progress")
library("tidyverse") 
library("dplyr")
library("TCGAbiolinks") 
library("stringr") 
library("TCGAutils") 
library('tibble') 
library('caret') 

data_dir <- "dataset"
RNA_sample_sheet_path <- paste0(data_dir, "/RNA_sample_sheet.tsv")

# Curated ARID1A interactors data from BioGRID
BioGRID_HSP_dir <- paste0(data_dir,'/BioGRID/BIOGRID-ORGANISM-Homo_sapiens.tab3.txt')

#BioGRID ARID1A PPI gene list extraction
biogrid <- read.delim(file = BioGRID_HSP_dir, sep = '\t')
biogrid <- biogrid %>% 
  filter(Official.Symbol.Interactor.A == "ARID1A" | Official.Symbol.Interactor.B == "ARID1A") %>% 
  filter(Experimental.System.Type == "genetic" | 
           Experimental.System %in% c('Two-hybrid','Affinity Capture-MS','Affinity Capture-RNA'))

PPI_list <- biogrid[,c('Official.Symbol.Interactor.A','Official.Symbol.Interactor.B','Experimental.System')]
PPI_list <- PPI_list %>% 
  mutate(Gene = ifelse(Official.Symbol.Interactor.A == 'ARID1A', Official.Symbol.Interactor.B, Official.Symbol.Interactor.A)) %>% 
  mutate(Gene = toupper(Gene)) %>% 
  select(c(Gene, Experimental.System))
PPI_list <- distinct(PPI_list, Gene, .keep_all = T)


#Pick UCEC mRNA data with PPI gene list
#sample number standard = mRNA count data
RNA_sample_sheet <- read.delim(RNA_sample_sheet_path, sep = "\t")
RNA_sample_sheet <- RNA_sample_sheet %>% 
  filter(Sample.Type == "Primary Tumor") %>% 
  filter(str_detect(Sample.ID, "-01A"))


#####################
## mRNA processing ##
#####################
PPI_frame <- data.frame(matrix(nrow = length(unique(RNA_sample_sheet[,"Case.ID"])), 
                                             ncol = length(PPI_list$Gene)))
rownames(PPI_frame) <- unique(RNA_sample_sheet[,"Case.ID"])
colnames(PPI_frame) <- PPI_list$Gene 
PPI_frame[] <- NaN

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
  for (PPI_component in PPI_list$Gene) {
    if (!(PPI_component %in% rna_file$gene_name)) {
      next
    }
    if (is.nan(PPI_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],PPI_component])) { ##
      PPI_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],PPI_component] <- rna_file[rna_file$gene_name == PPI_component,"unstranded"]
    }
    else {
      cur <- PPI_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],PPI_component]
      PPI_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],PPI_component] <- round((cur + rna_file[rna_file$gene_name == PPI_component,"unstranded"])/2)
    }
  }
}

PPI_frame <- PPI_frame[,!is.nan(colSums(PPI_frame))]

#feature scaling
#mRNA minmax scaling
scaled_PPI_frame <- data.frame(lapply(PPI_frame, function(x) scale(x, center = FALSE, scale = max(x, na.rm = TRUE)/1)))

#mRNA minmax scaling
process <- preProcess(PPI_frame, method = c("range"))
scaled_PPI_frame <- predict(process, PPI_frame)

#add protein
level4 <- read.table(paste0(data_dir,'/Protein/TCGA-UCEC-L4.csv'), sep = ',', header = T)
level4 <- level4 %>% 
  filter(str_detect(Sample_ID, "-01A")) %>%
  separate(Sample_ID,c('Case_ID','dump'), sep = '-01A', remove = FALSE, extra = 'merge')
level4 <- level4[,c("Case_ID","ARID1A")]
level4 <- rename(level4, "Row.names" = "Case_ID")
level4 <- level4[!(level4$Row.names=="TCGA-EY-A1GJ"),]
scaled_PPI_frame$Row.names <- rownames(scaled_PPI_frame)
scaled_PPI_frame <- merge(scaled_PPI_frame, level4, by = 'Row.names', all.x = T)

write.csv(scaled_PPI_frame, paste0(data_dir,'UCEC_ARID1A_BioGRID_PPI.csv'), sep = ',')
library("progress")
library("tidyverse") 
library("dplyr")
library("TCGAbiolinks")
library("stringr") 
library("TCGAutils") 
library('tibble') 

## The directory where you stored the queried dataset
## Save each data category as a subdirectory
## Sample sheet is "gdc_sample_sheet.YYYY-MM-DD.tsv" from GDC Data Portal
data_dir <- "dataset"
RNA_sample_sheet_path <- paste0(data_dir, "/RNA_sample_sheet.tsv")
CNV_sample_sheet_path <- paste0(data_dir, "/CNV_sample_sheet.tsv")
Met_sample_sheet_path <- paste0(data_dir, "/Met_sample_sheet.tsv")
miRNA_sample_sheet_path <- paste0(data_dir, "/miRNA_sample_sheet.tsv")
isomiR_sample_sheet_path <- paste0(data_dir, "/isomiR_sample_sheet.tsv")

#sample number baseline = mRNA count data
RNA_sample_sheet <- read.delim(RNA_sample_sheet_path, sep = "\t")
RNA_sample_sheet <- RNA_sample_sheet %>% 
  filter(Sample.Type == "Primary Tumor") %>% 
  filter(str_detect(Sample.ID, "-01A"))
Patient_data_frame <- data.frame(row.names = unique(RNA_sample_sheet[,"Case.ID"]))    

#####################
## mRNA processing ##
#####################
Patient_data_frame$RNA_count <- NaN
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
  if (is.nan(Patient_data_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],"RNA_count"])) { ##
    Patient_data_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],"RNA_count"] <- ARID1A_count
  }
  else {
    cur <- Patient_data_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],"RNA_count"]
    Patient_data_frame[RNA_sample_sheet[RNA_sample_sheet$File.Name == rna_file_name,"Case.ID"],"RNA_count"] <- round((cur + ARID1A_count)/2)
  }
}


####################
## CNV processing ##
####################
CNV_sample_sheet <- read.delim(CNV_sample_sheet_path, sep = "\t")
CNV_sample_sheet <- CNV_sample_sheet %>% 
  filter(str_detect(Sample.ID, "-01A"))

Patient_data_frame$CNV <- NaN
CNV_dir <- paste0(data_dir,"/CNV")
All_CNV_dir <- list.dirs(CNV_dir, recursive = FALSE)

pb <- progress_bar$new(total = length(All_CNV_dir)) #checking for loop progress
for (dir_name in All_CNV_dir) {
  pb$tick()
  cnv_file_name <- list.files(dir_name, pattern = "*.tsv")
  if (!(cnv_file_name %in% CNV_sample_sheet$File.Name)) {
    next
  }
  cnv_file <- read.delim(paste0(dir_name,"/",cnv_file_name),
                         sep = '\t',header = T,fill = T)
  ###CNV gain loss 
  ARID1A_CNV <- cnv_file[cnv_file$gene_name == "ARID1A", "copy_number"]
  Patient_data_frame[RNA_sample_sheet[CNV_sample_sheet$File.Name == cnv_file_name,"Case.ID"],"CNV"] <- ARID1A_CNV
}


############################
## Methylation processing ##
############################

#ARID1A associated CpG site list 
illumina <- read.csv(paste0(data_dir,"/illumina cpg probe id data.csv"),header = T, skip = 7)
illumina <- illumina[str_detect(illumina$UCSC_RefGene_Name,"ARID1A"),"IlmnID"]

##Add Beta Value of each site
Met_sample_sheet <- read.delim(Met_sample_sheet_path, sep = "\t")
Met_sample_sheet <- Met_sample_sheet %>% 
  filter(Sample.Type == "Primary Tumor") %>% 
  filter(str_detect(Sample.ID,"-01A"))
Patient_data_frame <- cbind(Patient_data_frame, setNames(lapply(illumina, function(x) x = NaN), illumina) )
Met_dir <- paste0(data_dir,"/Met")
All_Met_dir <- list.dirs(Met_dir, recursive = FALSE)

pb <- progress_bar$new(total = length(All_Met_dir)) #checking for loop progress
for (dir_name in All_Met_dir) {
  pb$tick()
  met_file_name <- list.files(dir_name, pattern = "*.level3betas.txt")
  if (!(met_file_name %in% Met_sample_sheet$File.Name)) {
    next
  }
  met_file <- read.table(paste0(dir_name,"/",met_file_name))
  met_file <- met_file[met_file$V1 %in% illumina,]
  rownames(met_file) <- met_file[,"V1"]
  Patient_data_frame[Met_sample_sheet[Met_sample_sheet$File.Name == met_file_name,"Case.ID"],rownames(met_file)] <- met_file[,"V2"]
}
## remove added useless row 
Patient_data_frame <- Patient_data_frame[rownames(Patient_data_frame) %in% unique(RNA_sample_sheet$Case.ID),]




################
## miRNA Data ##
################
miRNA_sample_sheet <- read.delim(miRNA_sample_sheet_path, sep = "\t")
miRNA_sample_sheet <- miRNA_sample_sheet %>% 
  filter(Sample.Type == "Primary Tumor") %>% 
  filter(str_detect(Sample.ID,"-01A"))
miRNA_dir <- paste0(data_dir,"/miRNA")
All_miRNA_dir <- list.dirs(miRNA_dir, recursive = FALSE)

isomiR_sample_sheet <- read.delim(isomiR_sample_sheet_path, sep = "\t")
isomiR_sample_sheet <- isomiR_sample_sheet %>% 
  filter(Sample.Type == "Primary Tumor") %>% 
  filter(str_detect(Sample.ID,"-01A"))
isomiR_dir <- paste0(data_dir,"/isomiR")
All_isomiR_dir <- list.dirs(isomiR_dir, recursive = FALSE)

miRDB_target <- read.delim(paste0(data_dir,"/miRDB_target.txt"),header = T, sep = '\t')
miRDB_target <- miRDB_target[,1:3]
arm_annotation_list <- read.delim(paste0(data_dir,"/miRNA_arm_feature_annotation.txt"),header = T, sep = '\t')
miRNA_dataframe <- data.frame(row.names = unique(RNA_sample_sheet[,"Case.ID"]))   
miRNA_dataframe[,miRDB_target$miRNA.Name] <- NaN

pb <- progress_bar$new(total = length(All_miRNA_dir))
for (miRNA in miRDB_target$miRNA.Name) {
  pb$tick()
  if (str_detect(miRNA, "-5p") | str_detect(miRNA,'-3p')) {
    miRNA_row <- which(arm_annotation_list == miRNA, arr.ind = T)[1]
    miRNA_col <- which(arm_annotation_list == miRNA, arr.ind = T)[2] + 1 #for MIMAT symbol column
    miRNA_symbol <- arm_annotation_list[miRNA_row, miRNA_col]
    for (dir_name in All_isomiR_dir) {
      isomiR_file_name <- list.files(dir_name, pattern = "*quantification.txt")
      isomiR_file <- read.table(paste0(dir_name,"/",isomiR_file_name))
      isomiR_rpm <- #top 3 isomiR to arm feature 
        sum(head(sort(as.numeric(isomiR_file[(which(isomiR_file == paste0("mature,",miRNA_symbol), arr.ind = T)[,1]),3]),decreasing = T),3))
      miRNA_dataframe[RNA_sample_sheet[isomiR_sample_sheet$File.Name == isomiR_file_name,"Case.ID"],miRNA] <- isomiR_rpm
    }
  }
  else {
    for (dir_name in All_miRNA_dir) {
      miRNA_file_name <- list.files(dir_name, pattern = "*quantification.txt")
      miRNA_file <- read.table(paste0(dir_name,"/",miRNA_file_name),header = T)
      miRNA_rpm <- miRNA_file[which(miRNA_file == tolower(miRNA), arr.ind = T)[1],"read_count"]
      miRNA_dataframe[RNA_sample_sheet[miRNA_sample_sheet$File.Name == miRNA_file_name,"Case.ID"],miRNA] <- miRNA_rpm
    }
  }
}

Patient_data_frame <- merge(x = Patient_data_frame,y = miRNA_dataframe, by = 'row.names', all.x = TRUE)
avg_rpm_1 <- data.frame(row.names = miRDB_target$miRNA.Name)
avg_rpm_1$avg <- 0
for (i in 1:length(miRDB_target$miRNA.Name)) {
  avg_rpm_1[i,"avg"] <- mean(Patient_data_frame[,miRDB_target$miRNA.Name[i]], na.rm = T)
}
avg_rpm_1$avg <- format(avg_rpm_1$avg, scientific = F, trim = T, digits = 2)


###################
## Mutation Data ##
###################

query <- GDCquery(
  project = "TCGA-UCEC", 
  data.category = "Simple Nucleotide Variation", 
  access = "open", 
  legacy = FALSE, 
  data.type = "Masked Somatic Mutation"
)
#GDCdownload(query)
SNV <- GDCprepare(query,directory = paste0(data_dir,'/Somatic'))
SNV <- SNV[SNV$Hugo_Symbol == "ARID1A",]

#filter features of interest 
filtered_SNV <- SNV[,c("case_id",'Start_Position','End_Position','Variant_Classification','HGVSp_Short','Consequence','IMPACT')]
pb <- progress_bar$new(total = length(filtered_SNV$case_id)) #checking for loop progress
for (UUID in filtered_SNV$case_id) {
  pb$tick()
  filtered_SNV[filtered_SNV$case_id == UUID,'case_id'] <- UUIDtoBarcode(id_vector = UUID)$submitter_id
}

#add features to Patient_data
filtered_SNV$Position <- paste(filtered_SNV$Start_Position, filtered_SNV$End_Position, filtered_SNV$Variant_Classification)
filtered_SNV <- filtered_SNV[,c('case_id','Position','Consequence','HGVSp_Short','IMPACT')]


##################
###File Export ###
##################
write.csv(Patient_data_frame, file = paste0(data_dir,"/TCGA_UCEC_summary.csv"))
write.csv(filtered_SNV, file = paste0(data_dir,"/TCGA_UCEC_mutation_summary.csv"))




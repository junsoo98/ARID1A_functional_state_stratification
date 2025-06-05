library("progress")
library("tidyverse") 
library("dplyr")
library("TCGAbiolinks") 
library("stringr") 
library("TCGAutils") 
library('tibble') 
library('caret') 

data_dir <- "dataset"
UCEC_summary <- read.csv(paste0(data_dir,"/TCGA_UCEC_summary.csv"), sep = ',')
UCEC_mutation <- read.csv(paste0(data_dir,"/TCGA_UCEC_mutation_summary.csv"), sep = ',')
UCEC_mutation[UCEC_mutation$Consequence == "missense_variant;splice_region_variant","Consequence"] <- "missense_variant"
UCEC_mutation[UCEC_mutation$Consequence == "stop_gained;splice_region_variant","Consequence"] <- "stop_gained"
rownames(UCEC_summary) <- UCEC_summary$Row.names


#mRNA minmax scaling
process <- preProcess(as.data.frame(UCEC_summary$RNA_count), method = c("range"))
scaled_RNA <- predict(process, as.data.frame(UCEC_summary$RNA_count))
UCEC_summary$scaled_RNA_count <- scaled_RNA$`UCEC_summary$RNA_count`


#CNV
UCEC_summary$CNV <- ifelse( UCEC_summary$CNV > 1, 0, 1)


#Methylation
UCEC_summary <- UCEC_summary[,-33] # remove ch.1.893623F because all NA
UCEC_summary[,5:32] <- mutate_all(UCEC_summary[,5:32], ~replace_na(.,0))

#miRNA normalizaiton
process <- preProcess(as.data.frame(UCEC_summary[,32:72]), method = c("range"))
UCEC_summary[,32:72] <- predict(process, as.data.frame(UCEC_summary[,32:72]))
miRNA_variance <- sapply(UCEC_summary[,32:72], var, na.rm = T)

#Mutation
mutation_features <- as.data.frame(table(UCEC_mutation$case_id)) #mutation number
row.names(mutation_features) <- mutation_features$Var1

mutation_features$FS_SG <- NaN #Frameshift, Stop gain
for (i in rownames(mutation_features)) {
  if ("stop_gained" %in% UCEC_mutation[(UCEC_mutation$case_id == i),"Consequence"] |
      "frameshift_variant" %in% UCEC_mutation[(UCEC_mutation$case_id == i),"Consequence"] |
      "missense_variant" %in% UCEC_mutation[(UCEC_mutation$case_id == i),"Consequence"]) {
    mutation_features[i,"FS_SG"] <- 1
  } else {
    mutation_features[i,"FS_SG"] <- 0.5
  }
}
mutation_features$IMPACT <- NaN #Mutation IMPACT score
for (i in rownames(mutation_features)) {
  if ("HIGH" %in% UCEC_mutation[(UCEC_mutation$case_id == i),"IMPACT"]) {
    mutation_features[i,"IMPACT"] <- 1
  } else if ("MODERATE" %in% UCEC_mutation[(UCEC_mutation$case_id == i),"IMPACT"]) {
    mutation_features[i,"IMPACT"] <- 0.5
  } else {
    mutation_features[i,"IMPACT"] <- 0
  }
}

mutation_features$Severe_Position <- NaN
UCEC_mutation$HGVSp_Short <- as.numeric(gsub(".*?([0-9]+).*", "\\1",  UCEC_mutation$HGVSp_Short))
#TCGA-AP-A0LE
for (i in rownames(mutation_features)) {
  mutation_position <- UCEC_mutation[(UCEC_mutation$case_id == i),"HGVSp_Short"]
  mutation_type <- UCEC_mutation[(UCEC_mutation$case_id == i),"Consequence"]
  mutation_type <- mutation_type %in% c("stop_gained","frameshift_variant")
  if (TRUE %in% mutation_type) {
    if (min(mutation_position[mutation_type]) < 1016) {
      mutation_features[i,"Severe_Position"] <- 1
    } else {mutation_features[i,"Severe_Position"] <- 0.75}
  } else {
    if (min(mutation_position) < 1016) {
      mutation_features[i,"Severe_Position"] <- 0.5
    } else {mutation_features[i,"Severe_Position"] <- 0}
  }
}

names(mutation_features)[names(mutation_features) == 'Var1'] <- "Row.names"
Final_UCEC <- merge(x = UCEC_summary, y = mutation_features, key = "Row.names", all.x = TRUE)
## remove miRNA features only with na 
to_remove <- c('hsa.miR.153.3p','hsa.miR.9.5p','hsa.miR.101.3p','hsa.miR.6511a.5p','hsa.miR.30c.5p','hsa.miR.548t.3p')
Final_UCEC <- Final_UCEC[, !(names(Final_UCEC) %in% to_remove)]

#Missing Vlaue imputation 
Final_UCEC[is.na(Final_UCEC$Freq),"Freq"] <- 0
Final_UCEC[is.na(Final_UCEC$FS_SG),"FS_SG"] <- 0
Final_UCEC[is.na(Final_UCEC$IMPACT),"IMPACT"] <- 0
Final_UCEC[is.na(Final_UCEC$Severe_Position),"Severe_Position"] <- 0

na_Final_UCEC <- which(is.na(Final_UCEC), arr.ind = T)
rownames(Final_UCEC) <- Final_UCEC$Row.names 

#Protein expression (TCPA)
level4 <- read.table(paste0(data_dir,'/Protein/TCGA-UCEC-L4.csv'), sep = ',', header = T)
level4 <- level4 %>% 
  filter(str_detect(Sample_ID, "-01A")) %>%
  separate(Sample_ID,c('Case_ID','dump'), sep = '-01A', remove = FALSE, extra = 'merge')
level4 <- level4[,c("Case_ID","ARID1A")]
level4 <- rename(level4, "Row.names" = "Case_ID")
level4 <- level4[!(level4$Row.names=="TCGA-EY-A1GJ"),]
Final_UCEC <- merge(Final_UCEC, level4, by = 'Row.names', all.x = T)


#Exporting processing 
rownames(Final_UCEC) <- Final_UCEC$Row.names
write.csv(Final_UCEC, paste0(data_dir,'TCGA_UCEC_scaled.csv'), sep = ',')

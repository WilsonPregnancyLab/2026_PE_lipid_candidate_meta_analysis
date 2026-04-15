# Packages (Run)
library(limma) #version 3.60.6
library(lumi) #version 2.56.0
library(stringr) #version 1.5.1
library(ggplot2) #version 3.5.1
library(gridExtra) #version 2.3
library(ggrepel) #version 0.9.5
library(minfi) #version 1.50.0
library(wateRmelon) #version 2.10.0
library(dplyr) #version 1.1.4

#Creating groups for sex stratification
#> colnames(affy_metadata)
#  [1] "X"                                   "title"                          
#  [3] "geo_accession"                       "apgar.score..1.min..ch1"        
#  [5] "apgar.score..5.min..ch1"             "attempted.vaginal.delivery.ch1" 
#  [7] "chorioamnionitis.diagnosis.ch1"      "molecule_ch1"                   
#  [9] "platform_id"                         "diagnosis.ch1"                  
# [11] "hellp.diagnosis.ch1"                 "infant.gender.ch1"              
# [13] "iugr.diagnosis.ch1"                  "maternal.age.ch1"               
# [15] "maternal.blood.type.ch1"             "maternal.bmi.ch1"               
# [17] "maternal.ethnicity.ch1"              "maximum.diastolic.bp.ch1"       
# [19] "maximum.systolic.bp.ch1"             "mean.umbilical.pi.ch1"          
# [21] "mean.uterine.pi.ch1"                 "mode.of.delivery.ch1"           
# [23] "mode.proteinuria.ch1"                "newborn.weight.z.score.ch1"     
# [25] "nicu.transfer.ch1"                   "placental.weight.z.score.ch1"   
# [27] "previous.hypertensive.pregnancy.ch1" "previous.miscarriage.ch1"       
# [29] "previous.nulliparity.ch1"            "umbilical.cord.diameter.ch1"    
# [31] "tissue.ch1"                          "source_name_ch1"                
# [33] "gestational_age"                     "file_path"  

full_affy_metadata <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/affymetrix_validation/metadata75010.csv")
affy_metadata <- full_affy_metadata[full_affy_metadata$file_path %in% pData(affymetrix_Lip_PE_final)$file_path, ]

as_factor <- function(col_names){affy_metadata[[col_names]] <<- as.factor(affy_metadata[[col_names]])}
col_list <- c("apgar.score..1.min..ch1", "chorioamnionitis.diagnosis.ch1", "attempted.vaginal.delivery.ch1", "chorioamnionitis.diagnosis.ch1", "diagnosis.ch1", "hellp.diagnosis.ch1", "infant.gender.ch1", "iugr.diagnosis.ch1", "maternal.age.ch1", "maternal.blood.type.ch1", "maternal.bmi.ch1", "maternal.ethnicity.ch1", "maximum.diastolic.bp.ch1", "maximum.systolic.bp.ch1", "mean.umbilical.pi.ch1", "mean.uterine.pi.ch1", "mode.of.delivery.ch1", "mode.proteinuria.ch1", "newborn.weight.z.score.ch1", "nicu.transfer.ch1", "placental.weight.z.score.ch1", "previous.hypertensive.pregnancy.ch1", "previous.miscarriage.ch1", "previous.nulliparity.ch1", "umbilical.cord.diameter.ch1", "gestational_age")
lapply(col_list, as_factor)

males <- affy_metadata[affy_metadata$infant.gender.ch1 == "M", ]
as_factor_males <- function(col_names){males[[col_names]] <<- as.factor(males[[col_names]])}
col_list <- c("apgar.score..1.min..ch1", "chorioamnionitis.diagnosis.ch1", "attempted.vaginal.delivery.ch1", "chorioamnionitis.diagnosis.ch1", "diagnosis.ch1", "hellp.diagnosis.ch1", "infant.gender.ch1", "iugr.diagnosis.ch1", "maternal.age.ch1", "maternal.blood.type.ch1", "maternal.bmi.ch1", "maternal.ethnicity.ch1", "maximum.diastolic.bp.ch1", "maximum.systolic.bp.ch1", "mean.umbilical.pi.ch1", "mean.uterine.pi.ch1", "mode.of.delivery.ch1", "mode.proteinuria.ch1", "newborn.weight.z.score.ch1", "nicu.transfer.ch1", "placental.weight.z.score.ch1", "previous.hypertensive.pregnancy.ch1", "previous.miscarriage.ch1", "previous.nulliparity.ch1", "umbilical.cord.diameter.ch1", "gestational_age")
lapply(col_list, as_factor_males)

females <- affy_metadata[affy_metadata$infant.gender.ch1 == "F", ]
as_factor_females <- function(col_names){females[[col_names]] <<- as.factor(females[[col_names]])}
col_list <- c("apgar.score..1.min..ch1", "chorioamnionitis.diagnosis.ch1", "attempted.vaginal.delivery.ch1", "chorioamnionitis.diagnosis.ch1", "diagnosis.ch1", "hellp.diagnosis.ch1", "infant.gender.ch1", "iugr.diagnosis.ch1", "maternal.age.ch1", "maternal.blood.type.ch1", "maternal.bmi.ch1", "maternal.ethnicity.ch1", "maximum.diastolic.bp.ch1", "maximum.systolic.bp.ch1", "mean.umbilical.pi.ch1", "mean.uterine.pi.ch1", "mode.of.delivery.ch1", "mode.proteinuria.ch1", "newborn.weight.z.score.ch1", "nicu.transfer.ch1", "placental.weight.z.score.ch1", "previous.hypertensive.pregnancy.ch1", "previous.miscarriage.ch1", "previous.nulliparity.ch1", "umbilical.cord.diameter.ch1", "gestational_age")
lapply(col_list, as_factor_females)


# Load in .csv file with significantly differentially expressed genes from RNA-seq analyses 

sig <- c("Decreased_RNA_Expression", "Increased_RNA_Expression", "Trending_Towards_Decreased_RNA_Expression", "Trending_Towards_Increased_RNA_Expression")

RNA_seq_autosomes_combined_sex <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/DESeq_rerun/RNA_Lip_DESeq_results_autosomes_combined_sex.csv") #dim 
RNA_seq_autosomes_combined_sex <- RNA_seq_autosomes_combined_sex[RNA_seq_autosomes_combined_sex$Expression_Status %in% sig, ] #dim 26, 12

RNA_seq_autosomes_F <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/DESeq_rerun/RNA_Lip_DESeq_results_autosomes_F.csv") #dim 
RNA_seq_autosomes_F <- RNA_seq_autosomes_F[RNA_seq_autosomes_F$Expression_Status %in% sig, ] #dim 1489, 12

RNA_seq_autosomes_M <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/DESeq_rerun/RNA_Lip_DESeq_results_autosomes_M.csv") #dim 
RNA_seq_autosomes_M <- RNA_seq_autosomes_M[RNA_seq_autosomes_M$Expression_Status %in% sig, ] #dim 1839, 12

RNA_seq_chrX_F <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/DESeq_rerun/RNA_Lip_DESeq_results_chrX_F.csv") #dim 
RNA_seq_chrX_F <- RNA_seq_chrX_F[RNA_seq_chrX_F$Expression_Status %in% sig, ] #dim 41, 12

RNA_seq_chrX_M <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/DESeq_rerun/RNA_Lip_DESeq_results_chrX_M.csv") #dim 
RNA_seq_chrX_M <- RNA_seq_chrX_M[RNA_seq_chrX_M$Expression_Status %in% sig, ] #dim 57, 12

RNA_seq_autosomes_combined_sex$ENS_ID <- sub("\\..*", "", RNA_seq_autosomes_combined_sex$Row.names)
RNA_seq_autosomes_F$ENS_ID <- sub("\\..*", "", RNA_seq_autosomes_F$Row.names)
RNA_seq_autosomes_M$ENS_ID <- sub("\\..*", "", RNA_seq_autosomes_M$Row.names)
RNA_seq_chrX_F$ENS_ID <- sub("\\..*", "", RNA_seq_chrX_F$Row.names)
RNA_seq_chrX_M$ENS_ID <- sub("\\..*", "", RNA_seq_chrX_M$Row.names)

#Probe annotations  
chrXprobes <- subset(fData(affymetrix_Lip_PE_final), fData(affymetrix_Lip_PE_final)$CHR == "X")
chrYprobes <- subset(fData(affymetrix_Lip_PE_final), fData(affymetrix_Lip_PE_final)$CHR == "Y")

# Read the filtered and normalized affymetrix signals, subset to only the significantly differentially expressed genes from the RNAseq analysis

affy_matrix <- exprs(affymetrix_Lip_PE_final)
affy_probe_data <- fData(affymetrix_Lip_PE_final)
affy_pheno <- pData(affymetrix_Lip_PE_final)

sig_probes_auto_combined_sex <- affy_probe_data[affy_probe_data$ENSEMBL %in% RNA_seq_autosomes_combined_sex$ENS_ID,] #dim 9, 6
sig_probes_auto_F <- affy_probe_data[affy_probe_data$ENSEMBL %in% RNA_seq_autosomes_F$ENS_ID,] #dim 1930, 6
sig_probes_auto_M <- affy_probe_data[affy_probe_data$ENSEMBL %in% RNA_seq_autosomes_M$ENS_ID,] #dim 1833, 6
sig_probes_chrX_F <- affy_probe_data[affy_probe_data$ENSEMBL %in% RNA_seq_chrX_F$ENS_ID,] #dim 54, 6
sig_probes_chrX_M <- affy_probe_data[affy_probe_data$ENSEMBL %in% RNA_seq_chrX_M$ENS_ID,] #dim 54, 6

affy_matrix_auto_combined_sex <- affy_matrix[rownames(affy_matrix) %in% sig_probes_auto_combined_sex$PROBEID,]
affy_matrix_auto_F <- affy_matrix[rownames(affy_matrix) %in% sig_probes_auto_F$PROBEID,]
affy_matrix_auto_M <- affy_matrix[rownames(affy_matrix) %in% sig_probes_auto_M$PROBEID,]
affy_matrix_chrX_F <- affy_matrix[rownames(affy_matrix) %in% sig_probes_chrX_F$PROBEID,]
affy_matrix_chrX_M <- affy_matrix[rownames(affy_matrix) %in% sig_probes_chrX_M$PROBEID,]


#Calculating the average Delta Beta values between Control and PE 
Controlmetadata_all <- affy_metadata[affy_metadata$diagnosis.ch1 == "non-PE", ]
Controlmetadata_F <- females[females$diagnosis.ch1 == "non-PE", ]
Controlmetadata_M <- males[males$diagnosis.ch1 == "non-PE", ]
# 72 control all (both males and females)
# 32 control females 
# 40 control males

PEmetadata_all <- affy_metadata[affy_metadata$diagnosis.ch1 == "PE", ]
PEmetadata_F <- females[females$diagnosis.ch1 == "PE", ]
PEmetadata_M <- males[males$diagnosis.ch1 == "PE", ]
# 74 PE all (both males and females)
# 37 PE females 
# 37 PE males

# Total Population Expression Table

#Autosomal Only Probes for Expression Table (these are not betas, just the signal)
expres_table <- function(matrix, PEmetadata, Controlmetadata){
    matrix_name <- deparse(substitute(matrix))
    suffix <- sub("affy_matrix_(.*)", "\\1", matrix_name)
    output_AvgExprs_name <- paste0("affy_AvgExpres_", suffix, ".csv")

    affy_PE <- as.data.frame(matrix[, PEmetadata$file_path])
    affy_CONT <- as.data.frame(matrix[, Controlmetadata$file_path])
    affy_CONT$AvgExprsCONT <- rowMeans(affy_CONT)
    affy_CONT$ProbeCONT <- rownames(affy_CONT)
    affy_PE$AvgExprsPE <- rowMeans(affy_PE)
    affy_PE$ProbePE <- rownames(affy_PE)
    #Merge All Table
    affy_avgExprs <- merge(affy_CONT[,c("AvgExprsCONT", "ProbeCONT")], affy_PE[,c("AvgExprsPE", "ProbePE")], by = "row.names")
    affy_avgExprs$deltaExprs <- affy_avgExprs$AvgExprsPE - affy_avgExprs$AvgExprsCONT
    rownames(affy_avgExprs) <- affy_avgExprs$ProbeCONT
    write.csv(affy_avgExprs, file = output_AvgExprs_name)
}

expres_table(affy_matrix_auto_combined_sex, PEmetadata_all, Controlmetadata_all)
expres_table(affy_matrix_auto_F, PEmetadata_F, Controlmetadata_F)
expres_table(affy_matrix_auto_M, PEmetadata_M, Controlmetadata_M)
expres_table(affy_matrix_chrX_F, PEmetadata_F, Controlmetadata_F)
expres_table(affy_matrix_chrX_M, PEmetadata_M, Controlmetadata_M)


#LINEAR MODELING

linear_modeling_combined_sex <- function (data, matrix, delta_file){
    data$gestational_age <- as.numeric(as.character(data$gestational_age))
    CONTvsPE_model <- model.matrix(~ diagnosis.ch1 + infant.gender.ch1 + gestational_age, data = data)
    CONTvsPE_fit <- lmFit(matrix, CONTvsPE_model)
    CONTvsPE_fit <- eBayes(CONTvsPE_fit)
    tt_CONTvsPE <- topTable(CONTvsPE_fit, n = Inf, adjust = "fdr", coef = "diagnosis.ch1PE")
    print(sum(tt_CONTvsPE$adj.P.Val < 0.05))
    tt_CONTvsPE$probes <- rownames(tt_CONTvsPE)
    tt_CONTvsPE$probes <- as.factor(tt_CONTvsPE$probes)
    delta_file <- read.csv(delta_file)
    delta_file$probes <- delta_file$ProbePE
    head(delta_file$probes)
    delta_file$probes <- as.factor(delta_file$probes)
    results <- merge(merge(tt_CONTvsPE, delta_file[, c("deltaExprs", "probes")], by = "probes"), 
            affy_probe_data[,c("PROBEID", "SYMBOL","GENENAME", "CHR", "ENSEMBL", "PMID")], by.x = "probes", by.y = "PROBEID"
    )
    sig_results <- results[results$adj.P.Val < 0.05, ]
    matrix_name <- deparse(substitute(matrix))
    suffix <- sub("affy_matrix_(.*)", "\\1", matrix_name)
    affy_results_name <- paste0("affy_results_", suffix, ".csv")
    tt_results_name <- paste0("tt_CONTvsPE_", suffix, ".csv")
    sig_results_name <- paste0("sig_affy_results_", suffix, ".csv")
    write.csv(results, file = affy_results_name)
    write.csv(tt_CONTvsPE, file = tt_results_name)
    write.csv(sig_results, file = sig_results_name)
} 
 
linear_modeling <- function (data, matrix, delta_file){
    data$gestational_age <- as.numeric(as.character(data$gestational_age))
    CONTvsPE_model <- model.matrix(~ diagnosis.ch1 + gestational_age, data = data)
    updated_matrix <- matrix[,colnames(matrix) %in% data$file_path]
    CONTvsPE_fit <- lmFit(updated_matrix, CONTvsPE_model)
    CONTvsPE_fit <- eBayes(CONTvsPE_fit)
    tt_CONTvsPE <- topTable(CONTvsPE_fit, n = Inf, adjust = "fdr", coef = "diagnosis.ch1PE")
    print(sum(tt_CONTvsPE$adj.P.Val < 0.05))
    tt_CONTvsPE$probes <- rownames(tt_CONTvsPE)
    tt_CONTvsPE$probes <- as.factor(tt_CONTvsPE$probes)
    delta_file <- read.csv(delta_file)
    delta_file$probes <- delta_file$ProbePE
    head(delta_file$probes)
    delta_file$probes <- as.factor(delta_file$probes)
    results <- merge(merge(tt_CONTvsPE, delta_file[, c("deltaExprs", "probes")], by = "probes"), 
            affy_probe_data[,c("PROBEID", "SYMBOL","GENENAME", "CHR", "ENSEMBL", "PMID")], by.x = "probes", by.y = "PROBEID"
    )
    sig_results <- results[results$adj.P.Val < 0.05, ]
    matrix_name <- deparse(substitute(matrix))
    suffix <- sub("affy_matrix_(.*)", "\\1", matrix_name)
    affy_results_name <- paste0("affy_results_", suffix, ".csv")
    tt_results_name <- paste0("tt_CONTvsPE_", suffix, ".csv")
    sig_results_name <- paste0("sig_affy_results_", suffix, ".csv")
    write.csv(results, file = affy_results_name)
    write.csv(tt_CONTvsPE, file = tt_results_name)
    write.csv(sig_results, file = sig_results_name)
} 


linear_modeling_combined_sex(affy_metadata, affy_matrix_auto_combined_sex, "affy_AvgExpres_auto_combined_sex.csv") #14 sig genes
linear_modeling(females, affy_matrix_auto_F, "affy_AvgExpres_auto_F.csv") #224 sig genes
linear_modeling(males, affy_matrix_auto_M, "affy_AvgExpres_auto_M.csv") #261 sig genes
linear_modeling(females, affy_matrix_chrX_F, "affy_AvgExpres_chrX_F.csv") #1 sig genes
linear_modeling(males, affy_matrix_chrX_M, "affy_AvgExpres_chrX_M.csv") #4 sig genes



#+ attempted.vaginal.delivery.ch1 + chorioamnionitis.diagnosis.ch1 + hellp.diagnosis.ch1 + iugr.diagnosis.ch1 + maximum.diastolic.bp.ch1 + maximum.systolic.bp.ch1 + mode.of.delivery.ch1 + newborn.weight.z.score.ch1 + nicu.transfer.ch1 + previous.miscarriage.ch1 + previous.nulliparity.ch1
# CONTvsPE_wholemodel_auto <- model.matrix(~ diagnosis.ch1 + infant.gender.ch1 + gestational_age , data = affy_metadata) 
# CONTvsPE_wholefit_auto <- lmFit(affy_matrix_auto_combined_sex, CONTvsPE_wholemodel_auto)
# CONTvsPE_wholefit_auto <- eBayes(CONTvsPE_wholefit_auto) #dim 96485, 22
# tt_CONTvsPE_whole_auto <- topTable(CONTvsPE_wholefit_auto, n = Inf, adjust = "fdr", coef = "diagnosis.ch1PE")
# print(sum(tt_CONTvsPE_whole_auto$adj.P.Val < 0.05)) #1
# tt_CONTvsPE_whole_auto$probes <- rownames(tt_CONTvsPE_whole_auto)
# affy_AvgExpres_auto_combined_sex <- read.csv("affy_AvgExpres_auto_combined_sex.csv")
# affy_AvgExpres_auto_combined_sex$probes <- row.names(affy_AvgExpres_auto_combined_sex) 
# tt_CONTvsPE_whole_auto$probes <- as.factor(tt_CONTvsPE_whole_auto$probes)
# affy_AvgExpres_auto_combined_sex$probes <- as.factor(affy_AvgExpres_auto_combined_sex$probes)
# affy_auto_combined_sex_results <- merge(merge(tt_CONTvsPE_whole_auto, affy_AvgExpres_auto_combined_sex[, c("deltaExprs","probes")], by = "probes"),
#                                affy_probe_data[,c("PROBEID", "SYMBOL","GENENAME", "CHR", "ENSEMBL", "PMID")], by.x = "probes", by.y = "PROBEID"
# )
# write.csv(tt_CONTvsPE_whole_auto, "./tt_CONTvsPE_whole_auto_study.csv")
# write.csv(placmet_wholepop_auto, "./placmet_wholepop_auto.csv") #129090 7
# auto_sig <- subset(placmet_wholepop_auto[placmet_wholepop_auto$adj.P.Val <0.05,]) #15 14 (8 unique probes but annotate to 15 different genes)
# write.csv(auto_sig, file = "./sig_all_autosomes_CONTvsPE_adjFunnorm.csv")

affy_auto_combined_sex <- read.csv ("affy_results_auto_combined_sex.csv")
affy_results_auto_F <- read.csv("affy_results_auto_F.csv")
affy_results_auto_M <- read.csv("affy_results_auto_M.csv")
affy_results_chrX_F <- read.csv("affy_results_chrX_F.csv")
affy_results_chrX_M <- read.csv("affy_results_chrX_M.csv")


#Whole data 
affy_auto_combined_sex$Expression_Status <- "Not_Biologically_Significant"
affy_auto_combined_sex$Expression_Status[affy_auto_combined_sex$logFC > 0.00 & affy_auto_combined_sex$adj.P.Val <0.05] <- "Trending Towards Increased Expression"
affy_auto_combined_sex$Expression_Status[affy_auto_combined_sex$logFC < 0.00 & affy_auto_combined_sex$adj.P.Val <0.05] <- "Trending Towards Decreased Expression"
affy_auto_combined_sex$Expression_Status[affy_auto_combined_sex$logFC > 1.00 & affy_auto_combined_sex$adj.P.Val <0.05] <- "Increased Expression"
affy_auto_combined_sex$Expression_Status[affy_auto_combined_sex$logFC < -1.00 & affy_auto_combined_sex$adj.P.Val <0.05] <- "Decreased Expression"

#Male data 
affy_results_auto_M$Expression_Status <- "Not_Biologically_Significant"
affy_results_auto_M$Expression_Status[affy_results_auto_M$logFC > 0.00 & affy_results_auto_M$adj.P.Val <0.05] <- "Trending Towards Increased Expression"
affy_results_auto_M$Expression_Status[affy_results_auto_M$logFC < 0.00 & affy_results_auto_M$adj.P.Val <0.05] <- "Trending Towards Decreased Expression"
affy_results_auto_M$Expression_Status[affy_results_auto_M$logFC > 1.00 & affy_results_auto_M$adj.P.Val <0.05] <- "Increased Expression"
affy_results_auto_M$Expression_Status[affy_results_auto_M$logFC < -1.00 & affy_results_auto_M$adj.P.Val <0.05] <- "Decreased Expression"

affy_results_chrX_M$Expression_Status <- "Not_Biologically_Significant"
affy_results_chrX_M$Expression_Status[affy_results_chrX_M$logFC > 0.00 & affy_results_chrX_M$adj.P.Val <0.05] <- "Trending Towards Increased Expression"
affy_results_chrX_M$Expression_Status[affy_results_chrX_M$logFC < 0.00 & affy_results_chrX_M$adj.P.Val <0.05] <- "Trending Towards Decreased Expression"
affy_results_chrX_M$Expression_Status[affy_results_chrX_M$logFC > 1.00 & affy_results_chrX_M$adj.P.Val <0.05] <- "Increased Expression"
affy_results_chrX_M$Expression_Status[affy_results_chrX_M$logFC < -1.00 & affy_results_chrX_M$adj.P.Val <0.05] <- "Decreased Expression"

#Female Data 
affy_results_auto_F$Expression_Status <- "Not_Biologically_Significant"
affy_results_auto_F$Expression_Status[affy_results_auto_F$logFC > 0.00 & affy_results_auto_F$adj.P.Val <0.05] <- "Trending Towards Increased Expression"
affy_results_auto_F$Expression_Status[affy_results_auto_F$logFC < 0.00 & affy_results_auto_F$adj.P.Val <0.05] <- "Trending Towards Decreased Expression"
affy_results_auto_F$Expression_Status[affy_results_auto_F$logFC > 1.00 & affy_results_auto_F$adj.P.Val <0.05] <- "Increased Expression"
affy_results_auto_F$Expression_Status[affy_results_auto_F$logFC < -1.00 & affy_results_auto_F$adj.P.Val <0.05] <- "Decreased Expression"

affy_results_chrX_F$Expression_Status <- "Not_Biologically_Significant"
affy_results_chrX_F$Expression_Status[affy_results_chrX_F$logFC > 0.00 & affy_results_chrX_F$adj.P.Val <0.05] <- "Trending Towards Increased Expression"
affy_results_chrX_F$Expression_Status[affy_results_chrX_F$logFC < 0.00 & affy_results_chrX_F$adj.P.Val <0.05] <- "Trending Towards Decreased Expression"
affy_results_chrX_F$Expression_Status[affy_results_chrX_F$logFC > 1.00 & affy_results_chrX_F$adj.P.Val <0.05] <- "Increased Expression"
affy_results_chrX_F$Expression_Status[affy_results_chrX_F$logFC < -1.00 & affy_results_chrX_F$adj.P.Val <0.05] <- "Decreased Expression"

combined_sex_auto_sig <- subset(affy_auto_combined_sex[affy_auto_combined_sex$adj.P.Val <0.05,])
female_auto_sig <- subset(affy_results_auto_F[affy_results_auto_F$adj.P.Val <0.05,]) #
male_auto_sig <- subset(affy_results_auto_M[affy_results_auto_M$adj.P.Val <0.05,]) #
female_chrX_sig <- subset(affy_results_chrX_F[affy_results_chrX_F$adj.P.Val <0.05,]) #
male_chrX_sig <- subset(affy_results_chrX_M[affy_results_chrX_M$adj.P.Val <0.05,]) #

#Volcano Plots "grey" (#no change in methylation), "#d02670"- (pink-Increased Expression), "#8a00c4"- (purple-Decreased Expression)

affy_auto_combined_sex$siglabel <- ifelse(affy_auto_combined_sex$Expression_Status %in% c("Increased Expression", "Decreased Expression"), affy_auto_combined_sex$SYMBOL, NA)
wholepop_auto <- ggplot(data = affy_auto_combined_sex, aes(x = logFC, y = -log10(adj.P.Val), col = Expression_Status)) + 
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 18)) +
  ylab("-log10(FDR)") +
  xlab("log2FoldChange") + 
  scale_y_continuous(breaks = seq(0, 12, by = 1), limits = c(0, 12)) +
  scale_x_continuous(breaks = seq(-2.5, 2.5, by = 0.5), limits = c(-2.5, 2.5)) +
  scale_color_manual(values = c("Decreased Expression" = "#8a00c4", "Increased Expression" = "#d02670", "Trending Towards Decreased Expression" = "grey", "Trending Towards Increased Expression" = "grey", "Not_Biologically_Significant" = "grey"),
                     guide = "none") +
  geom_point(shape = 19, alpha = 0.3, size = 3)+
  geom_vline(xintercept = c(-1, 1), col = "black", linetype = "dashed", linewidth = 0.75) +  
  geom_vline(xintercept = c(0), col = "black", linetype = "dashed", linewidth = 0.75) +
  geom_hline(yintercept = c(-log10(0.05)), col = "black", linetype = "dashed", linewidth = 0.75)  
  #geom_text_repel(aes(label=siglabel), na.rm = TRUE, max.overlaps = Inf, size = 4, segment.colour = 'grey50'  

  
affy_results_auto_M$siglabel <- ifelse(affy_results_auto_M$Expression_Status %in% c("Increased Expression", "Decreased Expression"), affy_results_auto_M$SYMBOL, NA)
male_auto <- ggplot(data = affy_results_auto_M, aes(x = logFC, y = -log10(adj.P.Val), col = Expression_Status)) + 
  theme_bw() +
  ylab(" ") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 18)) +
  xlab("log2FoldChange") +
  scale_y_continuous(breaks = seq(0, 12, by = 1), limits = c(0, 12)) +
  scale_x_continuous(breaks = seq(-2.5, 2.5, by = 0.5), limits = c(-2.5, 2.5)) +
  scale_color_manual(values = c("Decreased Expression" = "#8a00c4", "Increased Expression" = "#d02670", "Trending Towards Decreased Expression" = "grey", "Trending Towards Increased Expression" = "grey", "Not_Biologically_Significant" = "grey"),
                     guide = "none") + 
  geom_point(shape = 19, alpha = 0.4, size = 3) +
  geom_vline(xintercept = c(-1, 1), col = "black", linetype = "dashed", linewidth = 0.75) +
  geom_vline(xintercept = c(0), col = "black", linetype = "dashed", linewidth = 0.75) +
  geom_hline(yintercept = c(-log10(0.05)), col = "black", linetype = "dashed", linewidth = 0.75) +
  geom_text_repel(aes(label=siglabel), na.rm = TRUE, max.overlaps = Inf, size = 4, segment.colour = 'grey50')


affy_results_auto_F$siglabel <- ifelse(affy_results_auto_F$Expression_Status %in% c("Increased Expression", "Decreased Expression"), affy_results_auto_F$SYMBOL, NA)
female_auto <- ggplot(data = affy_results_auto_F, aes(x = logFC, y = -log10(adj.P.Val), col = Expression_Status)) + 
   theme_bw() +
  ylab("") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 18)) +
  xlab("log2FoldChange") +
  scale_y_continuous(breaks = seq(0, 12, by = 1), limits = c(0, 12)) +
  scale_x_continuous(breaks = seq(-2.5, 2.5, by = 0.5), limits = c(-2.5, 2.5)) +
  scale_color_manual(values = c("Decreased Expression" = "#8a00c4", "Increased Expression" = "#d02670", "Trending Towards Decreased Expression" = "grey", "Trending Towards Increased Expression" = "grey", "Not_Biologically_Significant" = "grey"),
                     guide = "none") +
  geom_point(shape = 19, alpha = 0.4, size = 3) +
  geom_vline(xintercept = c(-1, 1), col = "black", linetype = "dashed", linewidth = 0.75) +
  geom_vline(xintercept = c(0), col = "black", linetype = "dashed", linewidth = 0.75) +
  geom_hline(yintercept = c(-log10(0.05)), col = "black", linetype = "dashed", linewidth = 0.75) +
  geom_text_repel(aes(label=siglabel), na.rm = TRUE, max.overlaps = Inf, size = 4, segment.colour = 'grey50')



png("./all_autosome_vol_RNA_affy_panel.png", height = 9, width = 20, units = "in", res = 300)
grid.arrange(wholepop_auto, female_auto, male_auto, nrow = 1)
dev.off()

#plots of X chromosome 
affy_results_chrX_M$siglabel <- ifelse(affy_results_chrX_M$probes %in% male_chrX_sig$probes, affy_results_chrX_M$SYMBOL, NA)
male_X <- ggplot(data = affy_results_chrX_M, aes(x = logFC, y = -log10(adj.P.Val), col = Expression_Status)) +
  geom_point(shape = 19, alpha = 0.4, size = 3) +
  theme_bw() +
  ylab(" ") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 18)) +
  xlab("log2FoldChange") +
  scale_y_continuous(breaks = seq(0, 7, by = 0.5), limits = c(0, 7)) +
  scale_x_continuous(breaks = seq(-2.5, 2.5, by = 0.5), limits = c(-2.5, 2.5)) +
  scale_color_manual(values = c("Decreased Expression" = "#8a00c4", "Increased Expression" = "#d02670", "Trending Towards Decreased Expression" = "grey", "Trending Towards Increased Expression" = "grey", "Not_Biologically_Significant" = "grey"),
                     guide = "none") + 
geom_vline(xintercept = c(-1, 1), col = "black", linetype = "dashed", linewidth = 0.75) +
  geom_vline(xintercept = c(0), col = "black", linetype = "dashed", linewidth = 0.75) +
  geom_hline(yintercept = c(-log10(0.05)), col = "black", linetype = "dashed", linewidth = 0.75)


affy_results_chrX_F$siglabel <- ifelse(affy_results_chrX_F$probes %in% female_chrX_sig$probes, affy_results_chrX_F$SYMBOL, NA)
female_X <- ggplot(data = affy_results_chrX_F, aes(x = logFC, y = -log10(adj.P.Val), col = Expression_Status)) +
  geom_point(shape = 19, alpha = 0.4, size = 3) +
  theme_bw() +
  ylab("-log10(FDR)") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 18)) +
  xlab("log2FoldChange") +
  scale_y_continuous(breaks = seq(0, 7, by = 0.5), limits = c(0, 7)) +
  scale_x_continuous(breaks = seq(-2.5, 2.5, by = 0.5), limits = c(-2.5, 2.5)) +
  scale_color_manual(values = c("Decreased Expression" = "#8a00c4", "Increased Expression" = "#d02670", "Trending Towards Decreased Expression" = "grey", "Trending Towards Increased Expression" = "grey", "Not_Biologically_Significant" = "grey"),
                     guide = "none") + 
  #geom_text_repel(aes(label=siglabel), na.rm = TRUE, max.overlaps = Inf, size = 4, segment.colour = 'grey50'  + 
  geom_vline(xintercept = c(-1, 1), col = "black", linetype = "dashed", linewidth = 0.75) +
  geom_vline(xintercept = c(0), col = "black", linetype = "dashed", linewidth = 0.75) +
  geom_hline(yintercept = c(-log10(0.05)), col = "black", linetype = "dashed", linewidth = 0.75) 


png("./X_vol_RNA_affy_panel.png", height = 9, width = 15, units = "in", res = 300)
grid.arrange(female_X, male_X, nrow = 1)
dev.off()















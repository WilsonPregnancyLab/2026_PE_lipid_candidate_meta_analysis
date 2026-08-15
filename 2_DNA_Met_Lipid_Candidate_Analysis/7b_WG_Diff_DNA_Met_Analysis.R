# Packages (Run)
library(limma) #version 3.60.6
library(stringr) #version 1.5.1
library(ggplot2) #version 3.5.1
library(gridExtra) #version 2.3
library(ggrepel) #version 0.9.5
library(minfi) #version 1.50.0
library(wateRmelon) #version 2.10.0
library(dplyr) #version 1.1.4


# # Load in probe annotation files and lipid candidate gene file 

full_annotation_price38 <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_DNAm_Lipid_Candidate/Candidate_DNAm_Annotations/updated_full_annotation_region_vs_price38.csv") #dim 479561 10
full_annotation_price38$probes <- full_annotation_price38$CpG_ID
# # Subset annotation_price38 file so that it only includes my lipid genes
# lipid_candidate_probes_full <- annotation_price38[annotation_price38$gene %in% lipid_ensembl_list$gene_symbol, ] #113479 probes
annotation_price38 <- full_annotation_price38[!duplicated(full_annotation_price38$CpG_ID),] #414086 probes
lipid_candidate_probes <- read.csv("./lipid_candidate_probes.csv")

#Probe annotations 
chrXprobes <- subset(annotation_price38, annotation_price38$chr == "chrX") #10390 probes
chrYprobes <- subset(annotation_price38, annotation_price38$chr == "chrY") #380 probes

#Creating groups for sex stratification
metadata <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_DNAm_Lipid_Candidate/Candidate_Diff_DNAm_Analysis_2025/Metadata_Sheet_lipid_preeclampsia_excluded_removed.csv")
metadata$pathology_group <- as.factor(metadata$pathology_group)
metadata$Fetal_Sex <- as.factor(metadata$Fetal_Sex)
metadata$Sentrix_ID <- as.factor(metadata$Sentrix_ID)
metadata$Sentrix_Position <- as.factor(metadata$Sentrix_Position)
metadata$GSE_number <- as.factor(metadata$GSE_number)
metadata$geo_accession <- as.factor(metadata$geo_accession)
metadata$gestational_age <- as.numeric(as.character(metadata$gestational_age))

males <- subset(metadata, Fetal_Sex == "M")
males$pathology_group <- as.factor(males$pathology_group)
males$Fetal_Sex <- as.factor(males$Fetal_Sex)
males$Sentrix_ID <- as.factor(males$Sentrix_ID)
males$Sentrix_ID <- droplevels(males$Sentrix_ID)
males$Sentrix_Position <- as.factor(males$Sentrix_Position)
males$GSE_number <- as.factor(males$GSE_number)
males$gestational_age <- as.numeric(as.character(males$gestational_age))

females <- subset(metadata, Fetal_Sex == "F")
females$pathology_group <- as.factor(females$pathology_group)
females$Fetal_Sex <- as.factor(females$Fetal_Sex)
females$Sentrix_ID <- as.factor(females$Sentrix_ID)
females$Sentrix_ID <- droplevels(females$Sentrix_ID)
females$Sentrix_Position <- as.factor(females$Sentrix_Position)
females$GSE_number <- as.factor(females$GSE_number)
females$gestational_age <- as.numeric(as.character(females$gestational_age))

# Read the filtered and normalized RG sets
placmet_adjFunnorm_allfiltered <- readRDS("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_DNAm_Lipid_Candidate/Candidate_Diff_DNAm_Analysis_2025/placmet_adjFunnorm_allfiltered.rds") #dim 424750 172
placmet_adjFunnorm_filtbetas <- getBeta(placmet_adjFunnorm_allfiltered) #dim 424750 172
placmet_adjFunnorm_filtfun_F <- placmet_adjFunnorm_filtbetas[, females$Sample_Name]
# 83 females - gives all the GSMs in either Control or PE that are female 
placmet_adjFunnorm_filtfun_M <- placmet_adjFunnorm_filtbetas[, males$Sample_Name]
# 89 males - gives all the GSMs in either Control or PE that are male 

#Calculating the average Delta Beta values between Control and PE 
Controlmetadata_all <- subset(metadata, metadata$pathology_group == "Control") 
Controlmetadata_F <- subset(females, females$pathology_group == "Control")
Controlmetadata_M <- subset(males, males$pathology_group == "Control")
# 72 control all (both males and females)
# 35 control females 
# 37 control males
PEmetadata_all <- subset(metadata, metadata$pathology_group == "PE")
PEmetadata_F <- subset(females, females$pathology_group == "PE")
PEmetadata_M <- subset(males, males$pathology_group == "PE")
# 100 PE all
# 48 PE females 
# 52 PE males

# Total population Betas table 
# Autosomal Only Probes included for Beta Table

placmet_PE_all_autosomes <- as.data.frame(placmet_adjFunnorm_filtbetas[!rownames(placmet_adjFunnorm_filtbetas) %in% c(chrXprobes$CpG_ID, chrYprobes$CpG_ID), PEmetadata_all$Sample_Name]) #dim 415365 100
placmet_CONT_all_autosomes <- as.data.frame(placmet_adjFunnorm_filtbetas[!rownames(placmet_adjFunnorm_filtbetas) %in% c(chrXprobes$CpG_ID, chrYprobes$CpG_ID), Controlmetadata_all$Sample_Name]) #dim 415365 72
placmet_CONT_all_autosomes$AvgBetaCONT <- rowMeans(placmet_CONT_all_autosomes)
placmet_CONT_all_autosomes$ProbeCONT <- rownames(placmet_CONT_all_autosomes)
placmet_PE_all_autosomes$AvgBetaPE <- rowMeans(placmet_PE_all_autosomes)
placmet_PE_all_autosomes$ProbePE <- rownames(placmet_PE_all_autosomes)
# Merge All Table 
placmet_AllAvgbetas_autosomes <- merge(placmet_CONT_all_autosomes[,c("AvgBetaCONT","ProbeCONT")], placmet_PE_all_autosomes[,c("AvgBetaPE","ProbePE")], by = "row.names")
placmet_AllAvgbetas_autosomes$deltaB <- placmet_AllAvgbetas_autosomes$AvgBetaPE - placmet_AllAvgbetas_autosomes$AvgBetaCONT
rownames(placmet_AllAvgbetas_autosomes) <- placmet_AllAvgbetas_autosomes$ProbeCont
write.csv(placmet_AllAvgbetas_autosomes, "./placmet_AllAvgbetas_autosomes.csv")

# Male population Betas table 

#Autosomes
placmet_PE_M_autosomes <- as.data.frame(placmet_adjFunnorm_filtfun_M[!rownames(placmet_adjFunnorm_filtfun_M) %in% c(chrXprobes$CpG_ID,chrYprobes$CpG_ID), PEmetadata_M$Sample_Name])#dim: 415365, 52
placmet_CONT_M_autosomes <- as.data.frame(placmet_adjFunnorm_filtfun_M[!rownames(placmet_adjFunnorm_filtfun_M) %in% c(chrXprobes$CpG_ID,chrYprobes$CpG_ID), Controlmetadata_M$Sample_Name]) #dim: 415365, 37
# calculate Male average betas
placmet_CONT_M_autosomes$AvgBetaCONT <- rowMeans(placmet_CONT_M_autosomes)
placmet_CONT_M_autosomes$ProbeCONT <- rownames(placmet_CONT_M_autosomes)
placmet_PE_M_autosomes$AvgBetaPE <- rowMeans(placmet_PE_M_autosomes)
placmet_PE_M_autosomes$ProbePE <- rownames(placmet_PE_M_autosomes)
# Merge Male Table 
placmet_MaleAvgbetas_autosomes <- merge(placmet_CONT_M_autosomes[,c("AvgBetaCONT","ProbeCONT")], placmet_PE_M_autosomes[,c("AvgBetaPE","ProbePE")], by = "row.names")
placmet_MaleAvgbetas_autosomes$deltaB <- placmet_MaleAvgbetas_autosomes$AvgBetaPE - placmet_MaleAvgbetas_autosomes$AvgBetaCONT
rownames(placmet_MaleAvgbetas_autosomes) <- placmet_MaleAvgbetas_autosomes$ProbeCont
write.csv(placmet_MaleAvgbetas_autosomes, "./placmet_MaleAvgbetas_autosomes.csv")

#Male X Chromosome 
placmet_PE_M_X <-  as.data.frame(placmet_adjFunnorm_filtfun_M[rownames(placmet_adjFunnorm_filtfun_M) %in% chrXprobes$CpG_ID, PEmetadata_M$Sample_Name]) #dim: 9143, 52
placmet_CONT_M_X <- as.data.frame(placmet_adjFunnorm_filtfun_M[rownames(placmet_adjFunnorm_filtfun_M) %in% chrXprobes$CpG_ID, Controlmetadata_M$Sample_Name]) #dim: 9143, 37
# calculate Male X average betas
placmet_CONT_M_X$AvgBetaCONT <- rowMeans(placmet_CONT_M_X)
placmet_CONT_M_X$ProbeCONT <- rownames(placmet_CONT_M_X)
placmet_PE_M_X$AvgBetaPE <- rowMeans(placmet_PE_M_X)
placmet_PE_M_X$ProbePE <- rownames(placmet_PE_M_X)
# Merge Male X Table 
placmet_MaleAvgbetas_X <- merge(placmet_CONT_M_X[,c("AvgBetaCONT","ProbeCONT")], placmet_PE_M_X[,c("AvgBetaPE","ProbePE")], by = "row.names")
placmet_MaleAvgbetas_X$deltaB <- placmet_MaleAvgbetas_X$AvgBetaPE - placmet_MaleAvgbetas_X$AvgBetaCONT
rownames(placmet_MaleAvgbetas_X) <- placmet_MaleAvgbetas_X$ProbeCont
write.csv(placmet_MaleAvgbetas_X, "./placmet_MaleAvgbetas_X.csv")

#Male Y Chromosome (had 0 probes - didn't run, would copy previous code if had Male Y Chromosomes)
placmet_PE_M_Y <-  as.data.frame(placmet_adjFunnorm_filtfun_M[rownames(placmet_adjFunnorm_filtfun_M) %in% chrYprobes$CpG_ID, PEmetadata_M$Sample_Name]) #dim: 242, 52 
placmet_CONT_M_Y <- as.data.frame(placmet_adjFunnorm_filtfun_M[rownames(placmet_adjFunnorm_filtfun_M) %in% chrYprobes$CpG_ID, Controlmetadata_M$Sample_Name]) #dim: 242, 37
# calculate Male Y average betas
placmet_CONT_M_Y$AvgBetaCONT <- rowMeans(placmet_CONT_M_Y)
placmet_CONT_M_Y$ProbeCONT <- rownames(placmet_CONT_M_Y)
placmet_PE_M_Y$AvgBetaPE <- rowMeans(placmet_PE_M_Y)
placmet_PE_M_Y$ProbePE <- rownames(placmet_PE_M_Y)
# Merge Male Y Table 
placmet_MaleAvgbetas_Y <- merge(placmet_CONT_M_Y[,c("AvgBetaCONT","ProbeCONT")], placmet_PE_M_Y[,c("AvgBetaPE","ProbePE")], by = "row.names")
placmet_MaleAvgbetas_Y$deltaB <- placmet_MaleAvgbetas_Y$AvgBetaPE - placmet_MaleAvgbetas_Y$AvgBetaCONT
rownames(placmet_MaleAvgbetas_Y) <- placmet_MaleAvgbetas_Y$ProbeCont
write.csv(placmet_MaleAvgbetas_Y, "./placmet_MaleAvgbetas_Y.csv")


# Female population Betas table 

#Female Autosomes
placmet_PE_F_autosomes <-  as.data.frame(placmet_adjFunnorm_filtfun_F[!rownames(placmet_adjFunnorm_filtfun_F) %in% c(chrXprobes$CpG_ID,chrYprobes$CpG_ID), PEmetadata_F$Sample_Name]) #dim: 415365, 48
placmet_CONT_F_autosomes <- as.data.frame(placmet_adjFunnorm_filtfun_F[!rownames(placmet_adjFunnorm_filtfun_F) %in% c(chrXprobes$CpG_ID,chrYprobes$CpG_ID), Controlmetadata_F$Sample_Name]) #dim: 415365, 35
# calculate Female average betas
placmet_CONT_F_autosomes$AvgBetaCONT <- rowMeans(placmet_CONT_F_autosomes)
placmet_CONT_F_autosomes$ProbeCONT <- rownames(placmet_CONT_F_autosomes)
placmet_PE_F_autosomes$AvgBetaPE <- rowMeans(placmet_PE_F_autosomes)
placmet_PE_F_autosomes$ProbePE <- rownames(placmet_PE_F_autosomes)
# Merge Female Table 
placmet_FemaleAvgbetas_autosomes <- merge(placmet_CONT_F_autosomes[,c("AvgBetaCONT","ProbeCONT")], placmet_PE_F_autosomes[,c("AvgBetaPE","ProbePE")], by = "row.names")
placmet_FemaleAvgbetas_autosomes$deltaB <- placmet_FemaleAvgbetas_autosomes$AvgBetaPE - placmet_FemaleAvgbetas_autosomes$AvgBetaCONT
rownames(placmet_FemaleAvgbetas_autosomes) <- placmet_FemaleAvgbetas_autosomes$ProbeCont
write.csv(placmet_FemaleAvgbetas_autosomes, "./placmet_FemaleAvgbetas_autosomes.csv")

#Female X Chromosome 
placmet_PE_F_X <-  as.data.frame(placmet_adjFunnorm_filtfun_F[rownames(placmet_adjFunnorm_filtfun_F) %in% chrXprobes$CpG_ID, PEmetadata_F$Sample_Name])  #dim: 9143, 48
placmet_CONT_F_X <- as.data.frame(placmet_adjFunnorm_filtfun_F[rownames(placmet_adjFunnorm_filtfun_F) %in% chrXprobes$CpG_ID, Controlmetadata_F$Sample_Name])  #dim: 9143, 35
# calculate Female X average betas
placmet_CONT_F_X$AvgBetaCONT <- rowMeans(placmet_CONT_F_X)
placmet_CONT_F_X$ProbeCONT <- rownames(placmet_CONT_F_X)
placmet_PE_F_X$AvgBetaPE <- rowMeans(placmet_PE_F_X)
placmet_PE_F_X$ProbePE <- rownames(placmet_PE_F_X)
# Merge Female X Table 
placmet_FemaleAvgbetas_X <- merge(placmet_CONT_F_X[,c("AvgBetaCONT","ProbeCONT")], placmet_PE_F_X[,c("AvgBetaPE","ProbePE")], by = "row.names")
placmet_FemaleAvgbetas_X$deltaB <- placmet_FemaleAvgbetas_X$AvgBetaPE - placmet_FemaleAvgbetas_X$AvgBetaCONT
rownames(placmet_FemaleAvgbetas_X) <- placmet_FemaleAvgbetas_X$ProbeCont
write.csv(placmet_FemaleAvgbetas_X, "./placmet_FemaleAvgbetas_X.csv")

#m value inputs for linear models 
#Whole population 
placmet_adjFunnorm_MVal <- beta2m(placmet_adjFunnorm_filtbetas) #dim 424750, 172
write.csv(placmet_adjFunnorm_MVal, "./placmet_adjFunnorm_MVal.csv")
#whole population Mvalues - Autosomes
Mval_whole_auto <- placmet_adjFunnorm_MVal[!rownames(placmet_adjFunnorm_MVal) %in% c(chrXprobes$CpG_ID,chrYprobes$CpG_ID),] #424750, 172
write.csv(Mval_whole_auto, "./Mval_whole_auto.csv")

#Male Population Mvalues 
Mval_male <- beta2m(placmet_adjFunnorm_filtfun_M)
write.csv(Mval_male, "./MVal_male.csv")
#male autosomes
Mval_male_auto <- Mval_male[!rownames(Mval_male) %in% c(chrXprobes$CpG_ID,chrYprobes$CpG_ID),] #415365, 89
write.csv(Mval_male_auto, "./Mval_male_auto.csv")
#male XChr
Mval_male_X <- Mval_male[rownames(Mval_male) %in% chrXprobes$CpG_ID,] #9143, 89
write.csv(Mval_male_X, "./Mval_male_X.csv")
#male YChr
Mval_male_Y <- Mval_male[rownames(Mval_male) %in% chrYprobes$CpG_ID,] #242, 89
write.csv(Mval_male_Y, "./Mval_male_Y.csv")

#Female Population MVals 
Mval_female <- beta2m(placmet_adjFunnorm_filtfun_F)
write.csv(Mval_female, "./MVal_female.csv")
#female autosomes 
Mval_female_auto <- Mval_female[!rownames(Mval_female) %in% c(chrXprobes$CpG_ID,chrYprobes$CpG_ID),] #415365, 83
write.csv(Mval_female_auto, "./Mval_female_auto.csv")
#female XChr 
Mval_female_X <- Mval_female[rownames(Mval_female) %in% chrXprobes$CpG_ID,] #9143 83
write.csv(Mval_female_X, "./Mval_female_X.csv")

#LINEAR MODELING
#whole population autosomes

CONTvsPE_wholemodel_auto <- model.matrix(~ pathology_group + Fetal_Sex + GSE_number + gestational_age, data = metadata) 
CONTvsPE_wholefit_auto <- lmFit(Mval_whole_auto, CONTvsPE_wholemodel_auto)
CONTvsPE_wholefit_auto <- eBayes(CONTvsPE_wholefit_auto) 
tt_CONTvsPE_whole_auto <- topTable(CONTvsPE_wholefit_auto, n = Inf, adjust = "fdr", coef = "pathology_groupPE")
print(sum(tt_CONTvsPE_whole_auto$adj.P.Val < 0.05)) #0
tt_CONTvsPE_whole_auto$probes <- rownames(tt_CONTvsPE_whole_auto)
annotation_price38$probes <- annotation_price38$CpG_ID
names(placmet_AllAvgbetas_autosomes)[which(names(placmet_AllAvgbetas_autosomes) == "Row.names")] <- "probes"
tt_CONTvsPE_whole_auto$probes <- as.factor(tt_CONTvsPE_whole_auto$probes)
placmet_AllAvgbetas_autosomes$probes <- as.factor(placmet_AllAvgbetas_autosomes$probes)
placmet_wholepop_auto <- merge(merge(tt_CONTvsPE_whole_auto, placmet_AllAvgbetas_autosomes[, c("deltaB","probes")], by = "probes"),
                               full_annotation_price38[,c("gene", "probes", "chr","position", "region_overlap", "Closest_TSS_gene_name", "Closest_TSS_Transcript", "Closest_TSS_Pos", "Distance_Closest_TSS", "gene_id", "overlap")], by = "probes")
write.csv(tt_CONTvsPE_whole_auto, "./tt_CONTvsPE_whole_auto_study.csv")
write.csv(placmet_wholepop_auto, "./wg_placmet_wholepop_auto_rerun.csv") #129090 7
auto_sig <- subset(placmet_wholepop_auto[placmet_wholepop_auto$adj.P.Val <0.05,]) #15 14 (8 unique probes but annotate to 15 different genes)
write.csv(auto_sig, file = "./sig_all_autosomes_CONTvsPE_adjFunnorm.csv")

#whole population autosomes interaction term --> PF = pathology_group, fetal_sex

PF <- paste(metadata$pathology_group, metadata$Fetal_Sex, sep="_")
PF <- factor(PF, levels=c("Control_M","PE_M","Control_F","PE_F"))
CONTvsPE_wholemodel_auto <- model.matrix(~0 + PF + GSE_number + gestational_age, data = metadata) 
colnames(CONTvsPE_wholemodel_auto)[colnames(CONTvsPE_wholemodel_auto) == "PFControl_M"] <- "Control_M"
colnames(CONTvsPE_wholemodel_auto)[colnames(CONTvsPE_wholemodel_auto) == "PFPE_M"] <- "PE_M"
colnames(CONTvsPE_wholemodel_auto)[colnames(CONTvsPE_wholemodel_auto) == "PFControl_F"] <- "Control_F"
colnames(CONTvsPE_wholemodel_auto)[colnames(CONTvsPE_wholemodel_auto) == "PFPE_F"] <- "PE_F"
CONTvsPE_wholefit_auto <- lmFit(Mval_whole_auto, CONTvsPE_wholemodel_auto)

cont_matrix <- makeContrasts(
      PEinF=PE_F-Control_F,
      PEinM=PE_M-Control_M,
      Diff=(PE_F-Control_F)-(PE_M-Control_M),
      levels=CONTvsPE_wholefit_auto)
CONTvsPE_wholefit2_auto <- contrasts.fit(CONTvsPE_wholefit_auto, cont_matrix)
CONTvsPE_wholefit2_auto <- eBayes(CONTvsPE_wholefit2_auto)
tt_CONTvsPE_wholefit2_auto_F <- topTable(CONTvsPE_wholefit2_auto, n = Inf, adjust = "fdr", coef = "PEinF")
tt_CONTvsPE_wholefit2_auto_M <- topTable(CONTvsPE_wholefit2_auto, n = Inf, adjust = "fdr", coef = "PEinM")
tt_CONTvsPE_wholefit2_diff <- topTable(CONTvsPE_wholefit2_auto, n = Inf, adjust = "fdr", coef = "Diff")

print(sum(tt_CONTvsPE_wholefit2_auto_F$adj.P.Val < 0.05)) #0
print(sum(tt_CONTvsPE_wholefit2_auto_M$adj.P.Val < 0.05)) #3
print(sum(tt_CONTvsPE_wholefit2_diff$adj.P.Val < 0.05)) #0

tt_CONTvsPE_wholefit2_auto_F$probes <- tt_CONTvsPE_wholefit2_auto_F$X
annotation_price38$probes <- annotation_price38$CpG_ID
names(placmet_AllAvgbetas_autosomes)[which(names(placmet_AllAvgbetas_autosomes) == "Row.names")] <- "probes"
tt_CONTvsPE_wholefit2_auto_F$probes <- as.factor(tt_CONTvsPE_wholefit2_auto_F$probes)
placmet_AllAvgbetas_autosomes$probes <- as.factor(placmet_AllAvgbetas_autosomes$probes)
placmet_F_IT_auto <- merge(merge(tt_CONTvsPE_wholefit2_auto_F, placmet_AllAvgbetas_autosomes[, c("deltaB","probes")], by = "probes"),
                               annotation_price38[,c("gene", "probes", "chr","position", "region_overlap", "Closest_TSS_gene_name", "Closest_TSS_Transcript", "Closest_TSS_Pos", "Distance_Closest_TSS", "gene_id", "overlap")], by = "probes")
write.csv(tt_CONTvsPE_wholefit2_auto_F, "./tt_CONTvsPE_F_IT_auto_rerun.csv")
write.csv(placmet_F_IT_auto, "./wg_placmet_F_IT_auto_rerun.csv") #129090 15
auto_sig <- subset(placmet_F_IT_auto[placmet_F_IT_auto$adj.P.Val <0.05,]) #0

tt_CONTvsPE_wholefit2_auto_M$probes <- tt_CONTvsPE_wholefit2_auto_M$X
annotation_price38$probes <- annotation_price38$CpG_ID
names(placmet_AllAvgbetas_autosomes)[which(names(placmet_AllAvgbetas_autosomes) == "Row.names")] <- "probes"
tt_CONTvsPE_wholefit2_auto_M$probes <- as.factor(tt_CONTvsPE_wholefit2_auto_M$probes)
placmet_AllAvgbetas_autosomes$probes <- as.factor(placmet_AllAvgbetas_autosomes$probes)
placmet_M_IT_auto <- merge(merge(tt_CONTvsPE_wholefit2_auto_M, placmet_AllAvgbetas_autosomes[, c("deltaB","probes")], by = "probes"),
                               annotation_price38[,c("gene", "probes", "chr","position", "region_overlap", "Closest_TSS_gene_name", "Closest_TSS_Transcript", "Closest_TSS_Pos", "Distance_Closest_TSS", "gene_id", "overlap")], by = "probes")
write.csv(tt_CONTvsPE_wholefit2_auto_M, "./tt_CONTvsPE_M_IT_auto_rerun.csv")
write.csv(placmet_M_IT_auto, "./wg_placmet_M_IT_auto_rerun.csv") #129090 15
auto_sig <- subset(placmet_M_IT_auto[placmet_M_IT_auto$adj.P.Val <0.05,]) #0

#male autosomes
CONTvsPE_modelM <- model.matrix(~ pathology_group + GSE_number + gestational_age, data = males)
CONTvsPE_fitM_auto <- lmFit(Mval_male_auto, CONTvsPE_modelM) #warning: Partial NA coefficients for 1 probe(s)
CONTvsPE_fitM_auto <- eBayes(CONTvsPE_fitM_auto)
tt_CONTvsPE_M_all_auto <- topTable(CONTvsPE_fitM_auto, n = Inf, adjust = "fdr", coef = "pathology_groupPE")
print(sum(tt_CONTvsPE_M_all_auto$adj.P.Val < 0.05)) #19
tt_CONTvsPE_M_all_auto$probes <- rownames(tt_CONTvsPE_M_all_auto)
names(placmet_MaleAvgbetas_autosomes)[which(names(placmet_MaleAvgbetas_autosomes) == "Row.names")] <- "probes"
tt_CONTvsPE_M_all_auto$probes <- as.factor(tt_CONTvsPE_M_all_auto$probes)
placmet_MaleAvgbetas_autosomes$probes <- as.factor(placmet_MaleAvgbetas_autosomes$probes)
placmet_M_fulldata_auto <- merge(merge(tt_CONTvsPE_M_all_auto, placmet_MaleAvgbetas_autosomes[, c("deltaB","probes")], by = "probes"),
                               full_annotation_price38[,c("gene", "probes", "chr","position", "region_overlap", "Closest_TSS_gene_name", "Closest_TSS_Transcript", "Closest_TSS_Pos", "Distance_Closest_TSS", "gene_id", "overlap")], by = "probes")
write.csv(tt_CONTvsPE_M_all_auto, "./tt_CONTvsPE_M_all_auto_study.csv")
write.csv(placmet_M_fulldata_auto, "./wg_placmet_M_fulldata_auto_rerun.csv")
maleauto_sig <- subset(placmet_M_fulldata_auto[placmet_M_fulldata_auto$adj.P.Val <0.05,]) #0
write.csv(maleauto_sig, file = "./sig_Male_autosomes_CONTvsPE_adjFunnorm.csv")

#male X 
CONTvsPE_fitM_X <- lmFit(Mval_male_X, CONTvsPE_modelM) 
CONTvsPE_fitM_X <- eBayes(CONTvsPE_fitM_X)
tt_CONTvsPE_M_X <- topTable(CONTvsPE_fitM_X, n = Inf, adjust = "fdr", coef = "pathology_groupPE")
print(sum(tt_CONTvsPE_M_X$adj.P.Val < 0.05)) #0
tt_CONTvsPE_M_X$probes <- rownames(tt_CONTvsPE_M_X)
names(placmet_MaleAvgbetas_X)[which(names(placmet_MaleAvgbetas_X) == "Row.names")] <- "probes"
tt_CONTvsPE_M_X$probes <- as.factor(tt_CONTvsPE_M_X$probes)
placmet_MaleAvgbetas_X$probes <- as.factor(placmet_MaleAvgbetas_X$probes)
placmet_M_fulldata_X <- merge(merge(tt_CONTvsPE_M_X, placmet_MaleAvgbetas_X[, c("deltaB","probes")], by = "probes"),
                               full_annotation_price38[,c("gene", "probes", "chr","position", "region_overlap", "Closest_TSS_gene_name", "Closest_TSS_Transcript", "Closest_TSS_Pos", "Distance_Closest_TSS", "gene_id", "overlap")], by = "probes")
write.csv(tt_CONTvsPE_M_X, "./tt_CONTvsPE_M_X_study.csv")
write.csv(placmet_M_fulldata_X, "./wg_placmet_M_fulldata_X_rerun.csv")

#male Y
CONTvsPE_fitM_Y <- lmFit(Mval_male_Y, CONTvsPE_modelM) 
CONTvsPE_fitM_Y <- eBayes(CONTvsPE_fitM_Y)
tt_CONTvsPE_M_Y <- topTable(CONTvsPE_fitM_Y, n = Inf, adjust = "fdr", coef = "pathology_groupPE")
print(sum(tt_CONTvsPE_M_Y$adj.P.Val < 0.05)) #0
tt_CONTvsPE_M_Y$probes <- rownames(tt_CONTvsPE_M_Y)
names(placmet_MaleAvgbetas_Y)[which(names(placmet_MaleAvgbetas_Y) == "Row.names")] <- "probes"
tt_CONTvsPE_M_Y$probes <- as.factor(tt_CONTvsPE_M_Y$probes)
placmet_MaleAvgbetas_Y$probes <- as.factor(placmet_MaleAvgbetas_Y$probes)
placmet_M_fulldata_Y <- merge(merge(tt_CONTvsPE_M_Y, placmet_MaleAvgbetas_Y[, c("deltaB","probes")], by = "probes"),
                               full_annotation_price38[,c("gene", "probes", "chr","position", "region_overlap", "Closest_TSS_gene_name", "Closest_TSS_Transcript", "Closest_TSS_Pos", "Distance_Closest_TSS", "gene_id", "overlap")], by = "probes")
write.csv(tt_CONTvsPE_M_Y, "./tt_CONTvsPE_M_Y_study.csv")
write.csv(placmet_M_fulldata_Y, "./wg_placmet_M_fulldata_Y_rerun.csv")

#Female Autosomes 
CONTvsPE_modelF <- model.matrix(~ pathology_group + GSE_number + gestational_age, data = females)
CONTvsPE_fitF_auto <- lmFit(Mval_female_auto, CONTvsPE_modelF) 
CONTvsPE_fitF_auto <- eBayes(CONTvsPE_fitF_auto)
tt_CONTvsPE_F_all_auto <- topTable(CONTvsPE_fitF_auto, n = Inf, adjust = "fdr", coef = "pathology_groupPE")
print(sum(tt_CONTvsPE_F_all_auto$adj.P.Val < 0.05)) # 2
tt_CONTvsPE_F_all_auto$probes <- rownames(tt_CONTvsPE_F_all_auto)
names(placmet_FemaleAvgbetas_autosomes)[which(names(placmet_FemaleAvgbetas_autosomes) == "Row.names")] <- "probes"
tt_CONTvsPE_F_all_auto$probes <- as.factor(tt_CONTvsPE_F_all_auto$probes)
placmet_FemaleAvgbetas_autosomes$probes <- as.factor(placmet_FemaleAvgbetas_autosomes$probes)
placmet_F_fulldata_auto <- merge(merge(tt_CONTvsPE_F_all_auto, placmet_FemaleAvgbetas_autosomes[, c("deltaB","probes")], by = "probes"),
                               full_annotation_price38[,c("gene", "probes", "chr","position", "region_overlap", "Closest_TSS_gene_name", "Closest_TSS_Transcript", "Closest_TSS_Pos", "Distance_Closest_TSS", "gene_id", "overlap")], by = "probes")
write.csv(tt_CONTvsPE_F_all_auto, "./tt_CONTvsPE_F_all_auto_study.csv")
write.csv(placmet_F_fulldata_auto, "./wg_placmet_F_fulldata_auto_rerun.csv")
femaleauto_sig <- subset(placmet_F_fulldata_auto[placmet_F_fulldata_auto$adj.P.Val <0.05,])
write.csv(femaleauto_sig, file = "./sig_Female_autosomes_CONTvsPE_adjFunnorm.csv")

#Female X
CONTvsPE_fitF_X <- lmFit(Mval_female_X, CONTvsPE_modelF) 
CONTvsPE_fitF_X <- eBayes(CONTvsPE_fitF_X)
tt_CONTvsPE_F_X <- topTable(CONTvsPE_fitF_X, n = Inf, adjust = "fdr", coef = "pathology_groupPE")
print(sum(tt_CONTvsPE_F_X$adj.P.Val < 0.05)) #0
tt_CONTvsPE_F_X$probes <- rownames(tt_CONTvsPE_F_X)
names(placmet_FemaleAvgbetas_X)[which(names(placmet_FemaleAvgbetas_X) == "Row.names")] <- "probes"
tt_CONTvsPE_F_X$probes <- as.factor(tt_CONTvsPE_F_X$probes)
placmet_FemaleAvgbetas_X$probes <- as.factor(placmet_FemaleAvgbetas_X$probes)
placmet_F_fulldata_X <- merge(merge(tt_CONTvsPE_F_X, placmet_FemaleAvgbetas_X[, c("deltaB","probes")], by = "probes"),
                               full_annotation_price38[,c("gene", "probes", "chr","position", "region_overlap", "Closest_TSS_gene_name", "Closest_TSS_Transcript", "Closest_TSS_Pos", "Distance_Closest_TSS", "gene_id", "overlap")], by = "probes")
write.csv(tt_CONTvsPE_F_X, "./tt_CONTvsPE_F_X_study.csv")
write.csv(placmet_F_fulldata_X, "./wg_placmet_F_fulldata_X_rerun.csv")

#Plotting 
#Adding diff methylation information for colouring
#Read CSV files if plotting in a different session

placmet_wholepop_auto <- read.csv ("wg_placmet_wholepop_auto_rerun.csv")
placmet_M_fulldata_auto <- read.csv("wg_placmet_M_IT_auto_rerun.csv")
placmet_M_fulldata_X <- read.csv("wg_placmet_M_fulldata_X_rerun.csv")
placmet_F_fulldata_auto <- read.csv("wg_placmet_F_IT_auto_rerun.csv")
placmet_F_fulldata_X <- read.csv("wg_placmet_F_fulldata_X_rerun.csv")

#Whole data 
placmet_wholepop_auto$diffmethylation <- "Not_Biologically_Significant"
placmet_wholepop_auto$diffmethylation[placmet_wholepop_auto$deltaB > 0.00 & placmet_wholepop_auto$adj.P.Val <0.05] <- "Trending Towards Increased Methylation"
placmet_wholepop_auto$diffmethylation[placmet_wholepop_auto$deltaB < 0.00 & placmet_wholepop_auto$adj.P.Val <0.05] <- "Trending Towards Decreased Methylation"
placmet_wholepop_auto$diffmethylation[placmet_wholepop_auto$deltaB > 0.05 & placmet_wholepop_auto$adj.P.Val <0.05] <- "Increased Methylation"
placmet_wholepop_auto$diffmethylation[placmet_wholepop_auto$deltaB < -0.05 & placmet_wholepop_auto$adj.P.Val <0.05] <- "Decreased Methylation"
#Male data 
placmet_M_fulldata_auto$diffmethylation <- "Not_Biologically_Significant"
placmet_M_fulldata_auto$diffmethylation[placmet_M_fulldata_auto$deltaB > 0.00 & placmet_M_fulldata_auto$adj.P.Val <0.05] <- "Trending Towards Increased Methylation"
placmet_M_fulldata_auto$diffmethylation[placmet_M_fulldata_auto$deltaB < 0.00 & placmet_M_fulldata_auto$adj.P.Val <0.05] <- "Trending Towards Decreased Methylation"
placmet_M_fulldata_auto$diffmethylation[placmet_M_fulldata_auto$deltaB > 0.05 & placmet_M_fulldata_auto$adj.P.Val <0.05] <- "Increased Methylation"
placmet_M_fulldata_auto$diffmethylation[placmet_M_fulldata_auto$deltaB < -0.05 & placmet_M_fulldata_auto$adj.P.Val <0.05] <- "Decreased Methylation"
placmet_M_fulldata_X$diffmethylation <- "Not_Biologically_Significant"
placmet_M_fulldata_X$diffmethylation[placmet_M_fulldata_X$deltaB > 0.00 & placmet_M_fulldata_X$adj.P.Val <0.05] <- "Trending Towards Increased Methylation"
placmet_M_fulldata_X$diffmethylation[placmet_M_fulldata_X$deltaB < 0.00 & placmet_M_fulldata_X$adj.P.Val <0.05] <- "Trending Towards Decreased Methylation"
#Female Data 
placmet_F_fulldata_auto$diffmethylation <- "Not_Biologically_Significant"
placmet_F_fulldata_auto$diffmethylation[placmet_F_fulldata_auto$deltaB > 0.00 & placmet_F_fulldata_auto$adj.P.Val <0.05] <- "Trending Towards Increased Methylation"
placmet_F_fulldata_auto$diffmethylation[placmet_F_fulldata_auto$deltaB < 0.00 & placmet_F_fulldata_auto$adj.P.Val <0.05] <- "Trending Towards Decreased Methylation"
placmet_F_fulldata_auto$diffmethylation[placmet_F_fulldata_auto$deltaB > 0.05 & placmet_F_fulldata_auto$adj.P.Val <0.05] <- "Increased Methylation"
placmet_F_fulldata_auto$diffmethylation[placmet_F_fulldata_auto$deltaB < -0.05 & placmet_F_fulldata_auto$adj.P.Val <0.05] <- "Decreased Methylation"

placmet_F_fulldata_X$diffmethylation <- "Not_Biologically_Significant"
placmet_F_fulldata_X$diffmethylation[placmet_F_fulldata_X$deltaB > 0.00 & placmet_F_fulldata_X$adj.P.Val <0.05] <- "Trending Towards Increased Methylation"
placmet_F_fulldata_X$diffmethylation[placmet_F_fulldata_X$deltaB < 0.00 & placmet_F_fulldata_X$adj.P.Val <0.05] <- "Trending Towards Decreased Methylation"

# Comparison Table
library(dplyr) #version 1.1.4

# rename diffmethylation columns for better identification by group
# only ran the code for autosomes since there were no significant genes in the X chromosomes in fetal M and F populations
colnames(placmet_wholepop_auto)[colnames(placmet_wholepop_auto) == "diffmethylation"] <- "diffmethylation_whole"
colnames(placmet_F_fulldata_auto)[colnames(placmet_F_fulldata_auto) == "diffmethylation"] <- "diffmethylation_F"
colnames(placmet_M_fulldata_auto)[colnames(placmet_M_fulldata_auto) == "diffmethylation"] <- "diffmethylation_M"

# rename adj.P.Val columns for better identification by group
colnames(placmet_wholepop_auto)[colnames(placmet_wholepop_auto) == "adj.P.Val"] <- "adj.P.Val_whole"
colnames(placmet_F_fulldata_auto)[colnames(placmet_F_fulldata_auto) == "adj.P.Val"] <- "adj.P.Val_F"
colnames(placmet_M_fulldata_auto)[colnames(placmet_M_fulldata_auto) == "adj.P.Val"] <- "adj.P.Val_M"

# rename deltaB columns for better identification by group
colnames(placmet_wholepop_auto)[colnames(placmet_wholepop_auto) == "deltaB"] <- "deltaB_whole"
colnames(placmet_F_fulldata_auto)[colnames(placmet_F_fulldata_auto) == "deltaB"] <- "deltaB_F"
colnames(placmet_M_fulldata_auto)[colnames(placmet_M_fulldata_auto) == "deltaB"] <- "deltaB_M"

# subset groups to remove columns that are not biologically significant
sig_placmet_wholepop_auto <- subset(placmet_wholepop_auto[!placmet_wholepop_auto$diffmethylation_whole %in% c("Not_Biologically_Significant"), ])
sig_placmet_F_fulldata_auto <- subset(placmet_F_fulldata_auto[!placmet_F_fulldata_auto$diffmethylation_F %in% c("Not_Biologically_Significant"), ])
sig_placmet_M_fulldata_auto <- subset(placmet_M_fulldata_auto[!placmet_M_fulldata_auto$diffmethylation_M %in% c("Not_Biologically_Significant"), ])

# removed columns that were not being used for data analysis
sig_placmet_wholepop_auto_condensed <- sig_placmet_wholepop_auto[, (names(sig_placmet_wholepop_auto) %in% c("probes", "gene", "diffmethylation_whole", "adj.P.Val_whole", "deltaB_whole"))]
sig_placmet_F_fulldata_auto_condensed <- sig_placmet_F_fulldata_auto[, (names(sig_placmet_F_fulldata_auto) %in% c("probes", "gene", "diffmethylation_F", "adj.P.Val_F", "deltaB_F"))]
sig_placmet_M_fulldata_auto_condensed <- sig_placmet_M_fulldata_auto[, (names(sig_placmet_M_fulldata_auto) %in% c("probes", "Cgene", "diffmethylation_M", "adj.P.Val_M", "deltaB_M"))]

write.csv(sig_placmet_wholepop_auto_condensed, "sig_wholepop_auto_placmet.csv")
write.csv(sig_placmet_F_fulldata_auto_condensed, "sig_F_auto_placmet.csv")
write.csv(sig_placmet_M_fulldata_auto_condensed, "sig_M_auto_placmet.csv")

# One large table with all biologically significant genes in Female, Male and Whole Population Autosomes
all_sig_placmet <- merge(sig_placmet_wholepop_auto_condensed, sig_placmet_F_fulldata_auto_condensed, by = "probes", all = TRUE)
all_sig_placmet$gene <- coalesce(all_sig_placmet$gene.x, all_sig_placmet$gene.y)
all_sig_placmet <- all_sig_placmet[, !(names(all_sig_placmet) %in% c("gene.x", "gene.y"))]
write.csv(all_sig_placmet, "all_sig_placmet.csv")


wholeauto_sig <- subset(placmet_wholepop_auto[placmet_wholepop_auto$adj.P.Val_whole <0.05,]) #0
femaleauto_sig <- subset(placmet_F_fulldata_auto[placmet_F_fulldata_auto$adj.P.Val_F <0.05,]) #2
maleauto_sig <- subset(placmet_M_fulldata_auto[placmet_M_fulldata_auto$adj.P.Val_M <0.05,]) #16
wholeauto_bio_sig<- subset(placmet_wholepop_auto[placmet_wholepop_auto$adj.P.Val_whole <0.05 & (placmet_wholepop_auto$deltaB_whole < -0.05 | placmet_wholepop_auto$deltaB_whole > 0.05), ]) #0
femaleauto_bio_sig <- subset(placmet_F_fulldata_auto[placmet_F_fulldata_auto$adj.P.Val_F <0.05 & (placmet_F_fulldata_auto$deltaB_F < -0.05 | placmet_F_fulldata_auto$deltaB_F > 0.05), ]) #2
maleauto_bio_sig <- subset(placmet_M_fulldata_auto[placmet_M_fulldata_auto$adj.P.Val_M <0.05 & (placmet_M_fulldata_auto$deltaB_M < -0.05 | placmet_M_fulldata_auto$deltaB_M > 0.05), ]) #10



#Volcano Plots "grey" (#no change in methylation), "#d02670"- (pink-Increased Methylation), "#8a00c4"- (purple-Decreased Methylation)

placmet_wholepop_auto$siglabel <- ifelse(placmet_wholepop_auto$probes %in% wholeauto_bio_sig$probes, placmet_wholepop_auto$gene, NA)
wholepop_auto <- ggplot(data = placmet_wholepop_auto, aes(x = deltaB_whole, y = -log10(adj.P.Val_whole), col = diffmethylation_whole)) + 
  geom_point(shape = 19, alpha = 0.7, size = 5) + 
  theme_bw() +
  labs(title = "Combined-Sex Autosomal\nWhole-Genome Differential DNAm") +
  theme(plot.title = element_text(size = 17, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18)) +
  ylab("-log10(FDR)") +
  xlab("Delta Beta") + 
  scale_y_continuous(breaks = seq(0, 2.1, by = 0.5), limits = c(0, 2.1)) +
  scale_x_continuous(breaks = seq(-0.20, 0.20, by = 0.1), limits = c(-0.20, 0.20)) +
  scale_color_manual(values = c("grey","grey","grey","grey", "grey"),
                     guide = "none") +
  geom_vline(xintercept = c(-0.05,0.05), col = "black", linetype = "dashed", linewidth = 0.75) +
  geom_hline(yintercept = c(-log10(0.05)), col = "black", linetype = "dashed", linewidth = 0.75) 
 # geom_text_repel(aes(label=siglabel), na.rm = TRUE, max.overlaps = Inf, size = 4, segment.colour = 'grey50'  )
  
placmet_M_fulldata_auto$siglabel <- ifelse(placmet_M_fulldata_auto$probes %in% maleauto_bio_sig$probes, placmet_M_fulldata_auto$gene, NA)
male_auto <- ggplot(data = placmet_M_fulldata_auto, aes(x = deltaB_M, y = -log10(adj.P.Val_M), col = diffmethylation_M)) + 
geom_point(shape = 19, alpha = 0.7, size = 5) + 
  theme_bw() +
  labs(title = "Fetal Male Autosomal\nWhole-Genome Differential DNAm") +
  theme(plot.title = element_text(size = 17, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18)) +
          ylab(" ") +
  xlab("Delta Beta") +
  scale_y_continuous(breaks = seq(0, 2.1, by = 0.5), limits = c(0, 2.1)) +
  scale_x_continuous(breaks = seq(-0.20, 0.20, by = 0.1), limits = c(-0.20, 0.20)) +
  scale_color_manual(values = c("grey","grey","grey","grey", "grey"), 
                     guide = "none") + 
geom_vline(xintercept = c(-0.05,0.05), col = "black", linetype = "dashed", linewidth = 0.75) +
  geom_hline(yintercept = c(-log10(0.05)), col = "black", linetype = "dashed", linewidth = 0.75) 
#geom_text_repel(aes(label=siglabel), na.rm = TRUE, max.overlaps = Inf, size = 4, segment.colour = 'grey50')


placmet_F_fulldata_auto$siglabel <- ifelse(placmet_F_fulldata_auto$probes %in% femaleauto_bio_sig$probes, placmet_F_fulldata_auto$gene, NA)
female_auto <- ggplot(data = placmet_F_fulldata_auto, aes(x = deltaB_F, y = -log10(adj.P.Val_F), col = diffmethylation_F)) + 
geom_point(shape = 19, alpha = 0.7, size = 5) + 
  theme_bw() +
  labs(title = "Fetal Female Autosomal\nWhole-Genome Differential DNAm") +
  theme(plot.title = element_text(size = 17, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18)) +
        ylab("")+
  xlab("Delta Beta") +
  scale_y_continuous(breaks = seq(0, 2.1, by = 0.5), limits = c(0, 2.1)) +
  scale_x_continuous(breaks = seq(-0.20, 0.20, by = 0.1), limits = c(-0.20, 0.20)) +
  scale_color_manual(values = c("grey","grey","grey","grey", "grey"), 
                     guide = "none") +
  geom_vline(xintercept = c(-0.05,0.05), col = "black", linetype = "dashed", linewidth = 0.75) +
  geom_hline(yintercept = c(-log10(0.05)), col = "black", linetype = "dashed", linewidth = 0.75) 
#  geom_text_repel(aes(label=siglabel), na.rm = TRUE, max.overlaps = Inf, size = 4, segment.colour = 'grey50')


png("./wg_autosome_vol_adjFunnorm_panel_rerun.png", height = 9, width = 15, units = "in", res = 300)
grid.arrange(wholepop_auto, female_auto, male_auto, nrow = 1)
dev.off()

#plots of X chromosome 
male_X <- ggplot(data = placmet_M_fulldata_X, aes(x = deltaB, y = -log10(adj.P.Val), col = diffmethylation)) + 
  geom_vline(xintercept = c(-0.05,0.05), col = "black", linetype = "dashed", linewidth = 0.75) +
  geom_hline(yintercept = c(-log10(0.05)), col = "black", linetype = "dashed", linewidth = 0.75) +
  geom_point(shape = 19, alpha = 0.7, size = 5) + 
  theme_bw() +
  labs(title = "Fetal Female X-Chr\nWhole-Genome Differential DNAm") +
  theme(plot.title = element_text(size = 17, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18)) +
  ylab(" ") +
  xlab("Delta Beta") +
  scale_y_continuous(breaks = seq(0, 2.1, by = 0.5), limits = c(0, 2.1)) +
  scale_x_continuous(breaks = seq(-0.1, 0.1, by = 0.1), 
                     limits = c(-0.1,0.1))  +
  scale_color_manual(values = c("grey", "grey"), 
                     labels = c("Not Biologically Significant", "Increased Methylation"),
                     guide = "none")

female_X <- ggplot(data = placmet_F_fulldata_X, aes(x = deltaB, y = -log10(adj.P.Val), col = diffmethylation)) + 
  geom_vline(xintercept = c(-0.05,0.05), col = "black", linetype = "dashed", linewidth = 0.75) +
  geom_hline(yintercept = c(-log10(0.05)), col = "black", linetype = "dashed", linewidth = 0.75) +
  geom_point(shape = 19, alpha = 0.7, size = 5) + 
  theme_bw() +
  labs(title = "Fetal Male X-Chr\nWhole-Genome Differential DNAm") +
  theme(plot.title = element_text(size = 17, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18)) +
  ylab("-log10(FDR)") +
  xlab("Delta Beta") +
  scale_y_continuous(breaks = seq(0, 2.1, by = 0.5), limits = c(0, 2.1)) +
  scale_x_continuous(breaks = seq(-0.1, 0.1, by = 0.1), 
                     limits = c(-0.1,0.1)) +
  scale_color_manual(values = c("gray", "grey"), 
                     labels = c("Not Biologically Significant", "Increased Methylation"),
                     guide = "none")
png("./WG_X_vol_adjFunnorm_panel_rerun.png", height = 9, width = 15, units = "in", res = 300)
grid.arrange(female_X, male_X, nrow = 1)
dev.off()
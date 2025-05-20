# Packages (Run)
BiocManager::install("DMRcate")
library(DMRcate) 
library(missMethyl)
library(limma) #version 3.60.6

#Creating groups for sex stratification
metadata <- read.csv("/workspace/lab/wilsonslab/eyerk/2025_Lipid_Candidate_GSE_Info/GSE_metadata/Metadata_Sheet_lipid_preeclampsia_excluded_removed.csv")
metadata$pathology_group <- as.factor(metadata$pathology_group)
metadata$Fetal_Sex <- as.factor(metadata$Fetal_Sex)
metadata$Sentrix_ID <- as.factor(metadata$Sentrix_ID)
metadata$Sentrix_Position <- as.factor(metadata$Sentrix_Position)
metadata$GSE_number <- as.factor(metadata$GSE_number)
metadata$geo_accession <- as.factor(metadata$geo_accession)
metadata$gestational_age <- as.factor(metadata$gestational_age)

males <- subset(metadata, Fetal_Sex == "M")
males$pathology_group <- as.factor(males$pathology_group)
males$Fetal_Sex <- as.factor(males$Fetal_Sex)
males$Sentrix_ID <- as.factor(males$Sentrix_ID)
males$Sentrix_ID <- droplevels(males$Sentrix_ID)
males$Sentrix_Position <- as.factor(males$Sentrix_Position)
males$GSE_number <- as.factor(males$GSE_number)

females <- subset(metadata, Fetal_Sex == "F")
females$pathology_group <- as.factor(females$pathology_group)
females$Fetal_Sex <- as.factor(females$Fetal_Sex)
females$Sentrix_ID <- as.factor(females$Sentrix_ID)
females$Sentrix_ID <- droplevels(females$Sentrix_ID)
females$Sentrix_Position <- as.factor(females$Sentrix_Position)
females$GSE_number <- as.factor(females$GSE_number)

Controlmetadata_all <- subset(metadata, metadata$pathology_group == "Control") 
Controlmetadata_F <- subset(females, females$pathology_group == "Control")
Controlmetadata_M <- subset(males, males$pathology_group == "Control")

PEmetadata_all <- subset(metadata, metadata$pathology_group == "PE")
PEmetadata_F <- subset(females, females$pathology_group == "PE")
PEmetadata_M <- subset(males, males$pathology_group == "PE"

                       
# Read the filtered and normalized RG sets
placmet_adjFunnorm_allfiltered <- readRDS("/workspace/lab/wilsonslab/eyerk/2025_Thesis_Lipid_Candidate/R_entries/placmet_adjFunnorm_allfiltered.rds") #dim 329533 180
placmet_adjFunnorm_filtbetas_all <- getBeta(placmet_adjFunnorm_allfiltered)
placmet_adjFunnorm_filtbetas <- placmet_adjFunnorm_filtbetas_all[rownames(placmet_adjFunnorm_filtbetas_all) %in% lipid_candidate_probes$ID,] #dim 82401 180
placmet_adjFunnorm_filtfun_F <- placmet_adjFunnorm_filtbetas[, females$Sample_Name]
placmet_adjFunnorm_filtfun_M <- placmet_adjFunnorm_filtbetas[, males$Sample_Name]

# Combined Autosomes
combined_autosomes_beta <- placmet_adjFunnorm_filtbetas[!rownames(placmet_adjFunnorm_filtbetas) %in% c(chrXprobes$ID, chrYprobes$ID, NAprobes$ID), PEmetadata_all$Sample_Name] 
# Female Autosomes
F_autosomes_beta <- combined_autosomes[, females$Sample_Name]
F_X_beta <- placmet_adjFunnorm_filtbetas[rownames(placmet_adjFunnorm_filtbetas) %in% c(chrXprobes$ID, chrYprobes$ID, NAprobes$ID), PEmetadata_all$Sample_Name] 

# Male Autosomes
M_autosomes_beta <- combined_autosomes[, males$Sample_Name]
M_X_beta <- placmet_adjFunnorm_filtfun_M[rownames(placmet_adjFunnorm_filtfun_M) %in% chrXprobes$ID, PEmetadata_M$Sample_Name])

#DMR whole population autosomes

DMR_whole_auto_model <- model.matrix(~ pathology_group + Fetal_Sex + GSE_number + gestational_age, data = metadata) 
DMR_whole_auto_annot <- cpg.annotate("array", combined_autosomes_beta, what = "Beta", arraytype = "450K", analysis.type = "differential", design = DMR_whole_auto_model, fdr = 0.05, coef = 2)
DMR_whole_auto <- dmrcate(DMR_whole_auto_annot, lambda = 1000, C = 2)
results_range_whole_auto <- extractRanges(DMR_whole_auto, genome = "hg19")
wholeauto_bio_sig<- subset(results_range_whole_auto[results_range_whole_auto$HMFDR <0.05 & (results_range_whole_auto$meandiff < -0.05 | results_range_whole_auto$meandiff > 0.05), ]) #96


#GO Pathway Enrichment - MissMethyl
whole_auto_GO <- goregion(results_range_whole_auto[1:100], all.cpg = rownames(combined_autosomes_beta), collection = "GO", array.type = "450K")
whole_auto_GO <- whole_auto_GO[order(whole_auto_GO$P.DE),]























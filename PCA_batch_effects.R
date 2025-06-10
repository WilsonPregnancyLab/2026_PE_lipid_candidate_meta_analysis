# Packages (Run)
BiocManager::install("missMDA")
library(missMDA)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(magrittr)
library(wateRmelon)
library(ggplot2)
library(gridExtra) #version 2.3
library(ggrepel) #version 0.9.5

#Filtering unnormalized beta-values
placmet_noNorm <- readRDS("placmet_noNorm.RDS")

#Filtering bad probes
detp <- detectionP(placmet_noNorm)
probeInfo <- as.data.frame(cbind(IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Locations, IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Other, IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Manifest)) 
probeInfo$probeID <- rownames(probeInfo)
chrXprobes <- subset(probeInfo, probeInfo$chr == "chrX")
chrYprobes <- subset(probeInfo, probeInfo$chr == "chrY")
metadata <- read.csv("/workspace/lab/wilsonslab/eyerk/2025_Lipid_Candidate_GSE_Info/GSE_metadata/Metadata_Sheet_lipid_preeclampsia_excluded_removed.csv")

males <- subset(metadata, metadata$Fetal_Sex == "M")
females <- subset(metadata, metadata$Fetal_Sex == "F")
## For females, set detp in Y chromosomes to 0 (these probes do not bind to anything in female samples)
detp[rownames(detp) %in% chrYprobes$probeID, females$Sample_Name] <- 0
head(detp[rownames(detp) %in% chrYprobes$probeID, females$Sample_Name]) #Make sure all 0

# bad probes have detection p-value > 0.01
bad_detp <- detp > 0.01
#number of bad probes in >= 5% of samples
number_bad_detp <- print(sum(rowSums(bad_detp)>=(ncol(placmet_noNorm))*0.05)) #592
# missing betas >= 5% of samples 
avgbeta <- getBeta(placmet_noNorm)
bad_beta <- is.na(avgbeta)
number_bad_beta <- print(sum(rowSums(bad_beta)>=(ncol(placmet_noNorm))*0.05)) #17

# remove probes with bad p-values or missing betas >= 5% samples
badProbes <- rowSums(bad_detp)>=(ncol(placmet_noNorm))*0.05 | rowSums(bad_beta)>=(ncol(placmet_noNorm))*0.05 
placmet_noNorm_BPfilt <- avgbeta[!badProbes,]

# Remove SNP probes
# download full table from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL16304
price_anno <- read.delim("/workspace/lab/wilsonslab/eyerk/2025_Lipid_Candidate_GSE_Info/GSE_annotations/GPL16304-47833_no_legend.tsv",header=TRUE)
price_anno_SNP <- subset(price_anno, price_anno$n_target.CpG.SNP>0)
placmet_noNorm_SNPremoved <- placmet_noNorm_BPfilt[!rownames(placmet_noNorm_BPfilt) %in% price_anno_SNP$ID,]

# Remove Cross-Hybridizing probes
price_anno_Hyb <- subset(price_anno, price_anno$XY_Hits == "XY_YES" | price_anno$Autosomal_Hits == "A_YES")
placmet_noNorm_HybRemoved <- placmet_noNorm_SNPremoved[!rownames(placmet_noNorm_SNPremoved) %in% price_anno_Hyb$ID,]

# Remove Non-variable placental probes, #downloaded from https://github.com/redgar598/Tissue_Nonvariable_450K_CpGs
library(RCurl)
x <- getURL("https://raw.githubusercontent.com/redgar598/Tissue_Invariable_450K_CpGs/master/Invariant_Placenta_CpGs.csv")
write.csv (x,"Invariant_Placenta_CpGs.csv")
nonvarplac_anno <- read.csv("/workspace/lab/wilsonslab/lemairem/annotations/Invariant_Placenta_CpGs.csv", sep = ",") 
placmet_noNorm_nonvarremoved <- placmet_noNorm_HybRemoved[!rownames(placmet_noNorm_HybRemoved) %in% nonvarplac_anno$CpG,]
placmet_noNorm_allfiltered <- placmet_noNorm_nonvarremoved

saveRDS(placmet_noNorm_allfiltered, "placmet_noNorm_allfiltered.rds")

#Unnormalized PCA, have NA values so need to predict using missMDA
noNorm_pca_input <- t(placmet_noNorm_allfiltered)
estim <- estim_ncpPCA(noNorm_pca_input, scale = TRUE) #estimates number of components required for good model to predict NA values
estim$ncp #5, sub this value in the ncp = __ in the next line of code
sub_NA_noNorm <- imputePCA (noNorm_pca_input, ncp= 5, scale = TRUE) #subs NA values based on model
noNorm_pca <- prcomp(sub_NA_noNorm$completeObs, scale. = TRUE) #runs PCA

# Create a data frame with PC scores and metadata
noNorm_pca_df <- data.frame(noNorm_pca$x[, 1:2])  # PC1 and PC2
noNorm_pca_df$GSE_number <- metadata$GSE_number  
noNorm_pca_df$Fetal_Sex <- metadata$Fetal_Sex
noNorm_pca_df$pathology_group <- metadata$pathology_group

#Normalized PCA
placmet_adjFunnorm_allfiltered <- readRDS("/workspace/lab/wilsonslab/eyerk/2025_Thesis_Lipid_Candidate/R_entries/placmet_adjFunnorm_allfiltered.rds")
placmet_adjFunnorm_allfiltered_beta <- getBeta(placmet_adjFunnorm_allfiltered)
adjFunnorm_pca_input <- t(placmet_adjFunnorm_allfiltered_beta)
adjFunnorm_pca <- prcomp(adjFunnorm_pca_input, scale. = TRUE)

# Create a data frame with PC scores and metadata
adjFunnorm_pca_df <- data.frame(adjFunnorm_pca$x[, 1:2])  # PC1 and PC2
adjFunnorm_pca_df$GSE_number <- metadata$GSE_number
adjFunnorm_pca_df$Fetal_Sex <- metadata$Fetal_Sex
adjFunnorm_pca_df$pathology_group <- metadata$pathology_group

# Plot - GSE_number
PCA_noNorm <- ggplot(noNorm_pca_df, aes(x = PC1, y = PC2, color = GSE_number, fill = GSE_number)) +
  geom_point(size = 3) +
  theme_classic() +  
  scale_x_continuous(breaks = seq(-1000, 800, by = 250), limits = c(-1000, 800)) +
  scale_y_continuous(breaks = seq(-700, 600, by = 250), limits = c(-700, 600)) +
  geom_vline(xintercept = c(-0.05,0.05), col = "black", linetype = "dashed", linewidth = 0.75) +
  geom_hline(yintercept = c(-log10(0.05)), col = "black", linetype = "dashed", linewidth = 0.75) +
  ylab("PC2 (10.8%)") +
  xlab("PC1 (32.6%)") + 
  stat_ellipse(geom="polygon", level = 0.95, alpha = 0.2) +
  labs(title = "PCA Before Functional Normalization")

PCA_adjFunnorm <- ggplot(adjFunnorm_pca_df, aes(x = PC1, y = PC2, color = GSE_number, fill = GSE_number)) +
  geom_point(size = 3) +
  theme_classic() +
  scale_x_continuous(breaks = seq(-1000, 800, by = 250), limits = c(-1000, 800)) +
  scale_y_continuous(breaks = seq(-700, 600, by = 250), limits = c(-700, 600)) +  
  geom_vline(xintercept = c(-0.05,0.05), col = "black", linetype = "dashed", linewidth = 0.75) +
  geom_hline(yintercept = c(-log10(0.05)), col = "black", linetype = "dashed", linewidth = 0.75) +
  ylab("PC2 (9.7%)") +
  xlab("PC1 (12.2%)") + 
  stat_ellipse(geom="polygon", level = 0.95, alpha = 0.2) +
  labs(title = "PCA After Functional Normalization")

png("./PCA_Bef_Aft_Norm.png", height = 9, width = 30, units = "in", res = 300)
grid.arrange(PCA_noNorm, PCA_adjFunnorm, nrow = 1)
dev.off()

png("./PCA_noNorm.png", height = 9, width = 15, units = "in", res = 300)
grid.arrange(PCA_noNorm, nrow = 1)
dev.off()

png("./PCA_adjFunnorm.png", height = 9, width = 15, units = "in", res = 300)
grid.arrange(PCA_adjFunnorm, nrow = 1)
dev.off()

# Plot - Pathology Group 

PCA_noNorm_pg <- ggplot(noNorm_pca_df, aes(x = PC1, y = PC2, color = pathology_group, fill = pathology_group)) +
  geom_point(size = 3) +
  theme_classic() +  
  scale_x_continuous(breaks = seq(-1000, 800, by = 250), limits = c(-1000, 800)) +
  scale_y_continuous(breaks = seq(-700, 600, by = 250), limits = c(-700, 600)) +
  geom_vline(xintercept = c(-0.05,0.05), col = "black", linetype = "dashed", linewidth = 0.75) +
  geom_hline(yintercept = c(-log10(0.05)), col = "black", linetype = "dashed", linewidth = 0.75) +
  ylab("PC2 (10.8%)") +
  xlab("PC1 (32.6%)") + 
  stat_ellipse(geom="polygon", level = 0.95, alpha = 0.2) +
  labs(title = "PCA Before Functional Normalization")

PCA_adjFunnorm_pg <- ggplot(adjFunnorm_pca_df, aes(x = PC1, y = PC2, color = pathology_group, fill = pathology_group)) +
  geom_point(size = 3) +
  theme_classic() +
  scale_x_continuous(breaks = seq(-1000, 800, by = 250), limits = c(-1000, 800)) +
  scale_y_continuous(breaks = seq(-700, 600, by = 250), limits = c(-700, 600)) +  
  geom_vline(xintercept = c(-0.05,0.05), col = "black", linetype = "dashed", linewidth = 0.75) +
  geom_hline(yintercept = c(-log10(0.05)), col = "black", linetype = "dashed", linewidth = 0.75) +
  ylab("PC2 (9.7%)") +
  xlab("PC1 (12.2%)") + 
  stat_ellipse(geom="polygon", level = 0.95, alpha = 0.2) +
  labs(title = "PCA After Functional Normalization")

png("./PCA_Bef_Aft_Norm_pathology.png", height = 9, width = 30, units = "in", res = 300)
grid.arrange(PCA_noNorm_pg, PCA_adjFunnorm_pg, nrow = 1)
dev.off()

png("./PCA_noNorm_pathology.png", height = 9, width = 15, units = "in", res = 300)
grid.arrange(PCA_noNorm_pg, nrow = 1)
dev.off()

png("./PCA_adjFunnorm_pathology.png", height = 9, width = 15, units = "in", res = 300)
grid.arrange(PCA_adjFunnorm_pg, nrow = 1)
dev.off()


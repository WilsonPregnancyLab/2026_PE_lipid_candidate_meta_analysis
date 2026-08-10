# Packages (Run)
BiocManager::install("DMRcate")
library(DMRcate) #version 3.0.10
library(missMethyl) #version 1.38.0
library(limma) #version 3.60.6
library(ggplot2) #version 3.5.1
library(gridExtra) #version 2.3
library(ggrepel) #version 0.9.5

#Creating groups for sex stratification
metadata <- read.csv("/workspace/lab/wilsonslab/datalake-wilsonslab/PE_Lipid_Meta-analysis/2025_DNAm_Lipid_Candidate/Candidate_Diff_DNAm_Analysis_2025/Metadata_Sheet_lipid_preeclampsia_excluded_removed.csv")
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


# Subset annotation_price38 file so that it only includes my lipid genes
lipid_candidate_probes <- read.csv("/workspace/lab/wilsonslab/datalake-wilsonslab/PE_Lipid_Meta-analysis/2025_DNAm_Lipid_Candidate/Candidate_Diff_DNAm_Analysis_2025/lipid_candidate_probes.csv") #138524 17
chrXprobes <- subset(lipid_candidate_probes, lipid_candidate_probes$chr == "chrX") #2533 probes
chrYprobes <- subset(lipid_candidate_probes, lipid_candidate_probes$chr == "chrY") #0 probes
                       
# Read the filtered and normalized RG sets
placmet_adjFunnorm_allfiltered <- readRDS("/workspace/lab/wilsonslab/datalake-wilsonslab/PE_Lipid_Meta-analysis/2025_DNAm_Lipid_Candidate/Candidate_Diff_DNAm_Analysis_2025/placmet_adjFunnorm_allfiltered.rds") #dim 424750 172
placmet_adjFunnorm_filtbetas_all <- getBeta(placmet_adjFunnorm_allfiltered)
placmet_adjFunnorm_filtbetas <- placmet_adjFunnorm_filtbetas_all[rownames(placmet_adjFunnorm_filtbetas_all) %in% lipid_candidate_probes$CpG_ID,] #dim 124969 172
placmet_adjFunnorm_filtfun_F <- placmet_adjFunnorm_filtbetas[, females$Sample_Name] #124969 83
placmet_adjFunnorm_filtfun_M <- placmet_adjFunnorm_filtbetas[, males$Sample_Name] #124969 89

# Combined Autosomes
combined_auto_beta <- placmet_adjFunnorm_filtbetas[!rownames(placmet_adjFunnorm_filtbetas) %in% c(chrXprobes$CpG_ID, chrYprobes$CpG_ID), ] 
# Female Autosomes
F_auto_beta <- placmet_adjFunnorm_filtfun_F[!rownames(placmet_adjFunnorm_filtfun_F) %in% c(chrXprobes$CpG_ID, chrYprobes$CpG_ID), ] 
F_X_beta <- placmet_adjFunnorm_filtfun_F[rownames(placmet_adjFunnorm_filtfun_F) %in% c(chrXprobes$CpG_ID), ] 

# Male Autosomes
M_auto_beta <- placmet_adjFunnorm_filtfun_M[!rownames(placmet_adjFunnorm_filtfun_M) %in% c(chrXprobes$CpG_ID, chrYprobes$CpG_ID), ] 
M_X_beta <- placmet_adjFunnorm_filtfun_M[rownames(placmet_adjFunnorm_filtfun_M) %in% c(chrXprobes$CpG_ID), ]

#DMR whole population autosomes

DMR_whole_auto_model <- model.matrix(~ pathology_group + Fetal_Sex + GSE_number + gestational_age, data = metadata) 
DMR_whole_auto_annot <- cpg.annotate("array", combined_auto_beta, what = "Beta", arraytype = "450K", analysis.type = "differential", design = DMR_whole_auto_model, fdr = 0.05, coef = 2)
DMR_whole_auto <- dmrcate(DMR_whole_auto_annot, lambda = 1000, C = 2)
results_range_whole_auto <- extractRanges(DMR_whole_auto, genome = "hg19")
results_range_whole_auto$diffmethylation[results_range_whole_auto$meandiff > 0.00 & results_range_whole_auto$HMFDR <0.05] <- "Trending Towards Increased Methylation"
results_range_whole_auto$diffmethylation[results_range_whole_auto$meandiff < 0.00 & results_range_whole_auto$HMFDR <0.05] <- "Trending Towards Decreased Methylation"
results_range_whole_auto$diffmethylation[results_range_whole_auto$meandiff > 0.05 & results_range_whole_auto$HMFDR <0.05] <- "Increased Methylation"
results_range_whole_auto$diffmethylation[results_range_whole_auto$meandiff < -0.05 & results_range_whole_auto$HMFDR <0.05] <- "Decreased Methylation"
wholeauto_bio_sig<- subset(results_range_whole_auto[results_range_whole_auto$HMFDR <0.05 & (results_range_whole_auto$meandiff < -0.05 | results_range_whole_auto$meandiff > 0.05), ])
write.csv(results_range_whole_auto, "DMR_whole_auto_0.050.csv")
results_range_whole_auto_data <- as.data.frame(results_range_whole_auto)

#Plotting
#Volcano Plots
"grey" (#no change in methylation), "#d02670"- (pink-Increased Methylation), "#8a00c4"- (purple-Decreased Methylation)

whole_auto_DMR_plot <- ggplot(data = results_range_whole_auto_data, aes(x = meandiff, y = -log10(HMFDR), col = diffmethylation)) + 
  geom_vline(xintercept = c(-0.05,0.05), col = "black", linetype = "dashed", linewidth = 0.75) +
  geom_hline(yintercept = c(-log10(0.05)), col = "black", linetype = "dashed", linewidth = 0.75) +
  geom_point(shape = 19, alpha = 0.3, size = 3) + 
  theme_bw() +
  theme(axis.text = element_text(size = 12.5),
        axis.title = element_text(size = 14)) +
  ylab("-log10(FDR)") +
  xlab("Mean deltaB Across DMR") + 
  scale_y_continuous(breaks = seq(0, 5.5, by = 0.5), limits = c(0, 5.5)) +
  scale_x_continuous(breaks = seq(-0.10, 0.10, by = 0.05), limits = c(-0.10, 0.10)) +
  scale_color_manual(values = c("#8a00c4","#d02670","grey","grey", "grey"),
                     guide = "none")

  png("./wholepop_autos_vol_DMR.png", height = 9, width = 15, units = "in", res = 300)
  grid.arrange(whole_auto_DMR_plot, nrow = 1)
  dev.off()

#GO Pathway Enrichment - MissMethyl
whole_auto_GO <- goregion(results_range_whole_auto[1:100], all.cpg = rownames(combined_auto_beta), collection = "GO", array.type = "450K")
whole_auto_GO_sig<- subset(whole_auto_GO[whole_auto_GO$FDR <0.05,])
whole_auto_KEGG <- goregion(results_range_whole_auto[1:100], all.cpg = rownames(combined_auto_beta), collection = "KEGG", array.type = "450K")
whole_auto_KEGG_sig<- subset(whole_auto_KEGG[whole_auto_KEGG$FDR <0.05,])
write.csv(whole_auto_GO, "whole_auto_GO.csv")
write.csv(whole_auto_KEGG_sig, "whole_auto_KEGG_sig.csv")

#Function version
  
DMR_Analysis <- function(beta, data){
  print("DMR_Analysis")
  if (identical(data, metadata)) {
    model <- model.matrix(~ pathology_group + Fetal_Sex + GSE_number + gestational_age, data = data)
  } else {
    model <- model.matrix(~ pathology_group + GSE_number + gestational_age, data = data)
  }  
  model_annot <- cpg.annotate("array", beta, what = "Beta", arraytype = "450K", analysis.type = "differential", design = model, fdr = 0.05, coef = "pathology_groupPE")

  dmr_model_auto <- dmrcate(model_annot, lambda=1000, C=2)
  results_ranges <- extractRanges(dmr_model_auto, genome = "hg19")
  results_ranges$diffmethylation[results_ranges$meandiff > 0.00 & results_ranges$HMFDR <0.05] <- "Trending Towards Increased Methylation"
  results_ranges$diffmethylation[results_ranges$meandiff < 0.00 & results_ranges$HMFDR <0.05] <- "Trending Towards Decreased Methylation"
  results_ranges$diffmethylation[results_ranges$meandiff > 0.05 & results_ranges$HMFDR <0.05] <- "Increased Methylation"
  results_ranges$diffmethylation[results_ranges$meandiff < -0.05 & results_ranges$HMFDR <0.05] <- "Decreased Methylation"
  Bval_name <- deparse(substitute(beta))
  filename <- paste0(Bval_name, "_results_ranges.csv")
  write.csv(results_ranges, file = filename)
  results_ranges_bio_sig<- subset(results_ranges[results_ranges$HMFDR <0.05 & (results_ranges$meandiff < -0.05 | results_ranges$meandiff > 0.05), ])
  filename_bio_sig <- paste0(Bval_name, "_results_ranges_bio_sig.csv")
  write.csv(results_ranges_bio_sig, file = filename_bio_sig)

  print("DMR_Volcano_Plot")
  results_range_data <- as.data.frame(results_ranges)
  DMR_plot <- ggplot(data = results_range_data, aes(x = meandiff, y = -log10(HMFDR), col = diffmethylation)) + 
    geom_vline(xintercept = c(-0.05,0.05), col = "black", linetype = "dashed", linewidth = 0.75) +
    geom_hline(yintercept = c(-log10(0.05)), col = "black", linetype = "dashed", linewidth = 0.75) +
    geom_point(shape = 19, alpha = 0.3, size = 3) + 
    theme_bw() +
    theme(axis.text = element_text(size = 12.5),
          axis.title = element_text(size = 14)) +
    ylab("-log10(FDR)") +
    xlab("DMR Mean DeltaB") + 
    #scale_y_continuous(breaks = seq(0, 5.5, by = 0.5), limits = c(0, 5.5)) +
    #scale_x_continuous(breaks = seq(-0.10, 0.10, by = 0.05), limits = c(-0.10, 0.10)) +
    scale_color_manual(values = c("#8a00c4","#d02670","grey","grey", "grey"),
                       guide = "none")

  filename_vol_plot <- paste0(Bval_name, "_vol_DMR.png")
  png(file = filename_vol_plot, height = 9, width = 15, units = "in", res = 300)
  grid.arrange(DMR_plot, nrow = 1)
  dev.off()

  print("DMR_GO_ENRICHMENT")
  sig_regions <- results_ranges[results_ranges$HMFDR < 0.05]
  enrichment_GO <- goregion(sig_regions, all.cpg = rownames(beta), collection = "GO", array.type = "450K")
  enrichment_GO <- enrichment_GO[order(enrichment_GO$P.DE),]
  filename_GO <- paste0(Bval_name, "_DMR_enrichment_GO.csv")
  write.csv(enrichment_GO, file = filename_GO)}



whole_auto_DMR <- DMR_Analysis(combined_auto_beta, metadata)
#5 ind sig probes, no genes annotated to sig CpGs

male_auto_DMR <- DMR_Analysis(M_auto_beta, males)
#no ind sig probes, no genes annotated to sig CpGs

male_X_DMR <- DMR_Analysis(M_X_beta, males)
#no ind sig probes, no genes annotated to sig CpGs

female_auto_DMR <- DMR_Analysis(F_auto_beta, females)
#1 ind sig probes, 1 gene annotated to sig CpGs

female_X_DMR <- DMR_Analysis(F_X_beta, females)
#no ind sig probes, no genes annotated to sig CpGs














## 11b - UpSet plot - overlaps in DMRs between all 6 models
# get DMR names (chr_start_end format) for each model

# Packages (Run)
BiocManager::install(c("UpSetR", "tidyverse", "RColorBrewer", "gridExtra", "ggplot2"))
library(UpSetR) #version 1.4.0
library(tidyverse) #version 2.0.0
library(RColorBrewer) #version 1.1.3
library(gridExtra) #version 2.3

DNAm_Lip_auto_F <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_DNAm_Lipid_Candidate/Candidate_Diff_DNAm_Analysis_2025/placmet_F_fulldata_auto_rerun.csv")
DNAm_Lip_auto_M <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_DNAm_Lipid_Candidate/Candidate_Diff_DNAm_Analysis_2025/placmet_M_fulldata_auto_rerun.csv")
#Male data 
DNAm_Lip_auto_M$diffmethylation <- "Not_Biologically_Significant"
DNAm_Lip_auto_M$diffmethylation[DNAm_Lip_auto_M$deltaB > 0.00 & DNAm_Lip_auto_M$adj.P.Val <0.05] <- "Trending Towards Increased Methylation"
DNAm_Lip_auto_M$diffmethylation[DNAm_Lip_auto_M$deltaB < 0.00 & DNAm_Lip_auto_M$adj.P.Val <0.05] <- "Trending Towards Decreased Methylation"
DNAm_Lip_auto_M$diffmethylation[DNAm_Lip_auto_M$deltaB > 0.05 & DNAm_Lip_auto_M$adj.P.Val <0.05] <- "Increased Methylation"
DNAm_Lip_auto_M$diffmethylation[DNAm_Lip_auto_M$deltaB < -0.05 & DNAm_Lip_auto_M$adj.P.Val <0.05] <- "Decreased Methylation"
#Female Data 
DNAm_Lip_auto_F$diffmethylation <- "Not_Biologically_Significant"
DNAm_Lip_auto_F$diffmethylation[DNAm_Lip_auto_F$deltaB > 0.00 & DNAm_Lip_auto_F$adj.P.Val <0.05] <- "Trending Towards Increased Methylation"
DNAm_Lip_auto_F$diffmethylation[DNAm_Lip_auto_F$deltaB < 0.00 & DNAm_Lip_auto_F$adj.P.Val <0.05] <- "Trending Towards Decreased Methylation"
DNAm_Lip_auto_F$diffmethylation[DNAm_Lip_auto_F$deltaB > 0.05 & DNAm_Lip_auto_F$adj.P.Val <0.05] <- "Increased Methylation"
DNAm_Lip_auto_F$diffmethylation[DNAm_Lip_auto_F$deltaB < -0.05 & DNAm_Lip_auto_F$adj.P.Val <0.05] <- "Decreased Methylation"

RNA_Lip_auto_F <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/DESeq_rerun/RNA_Lip_DESeq_results_autosomes_F.csv")
RNA_Lip_auto_M <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/DESeq_rerun/RNA_Lip_DESeq_results_autosomes_M.csv")

affy_Lip_auto_F <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/affymetrix_validation_rerun/affy_results_auto_F.csv")
affy_Lip_auto_M <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/affymetrix_validation_rerun/affy_results_auto_M.csv")
#Add Significance Label - Male Data
affy_Lip_auto_M$Expression_Status <- "Not_Biologically_Significant"
affy_Lip_auto_M$Expression_Status[affy_Lip_auto_M$logFC > 0.00 & affy_Lip_auto_M$adj.P.Val <0.05] <- "Trending Towards Increased Expression"
affy_Lip_auto_M$Expression_Status[affy_Lip_auto_M$logFC < 0.00 & affy_Lip_auto_M$adj.P.Val <0.05] <- "Trending Towards Decreased Expression"
affy_Lip_auto_M$Expression_Status[affy_Lip_auto_M$logFC > 0.05 & affy_Lip_auto_M$adj.P.Val <0.05] <- "Increased Expression"
affy_Lip_auto_M$Expression_Status[affy_Lip_auto_M$logFC < -0.05 & affy_Lip_auto_M$adj.P.Val <0.05] <- "Decreased Expression"
#Female Data 
affy_Lip_auto_F$Expression_Status <- "Not_Biologically_Significant"
affy_Lip_auto_F$Expression_Status[affy_Lip_auto_F$logFC > 0.00 & affy_Lip_auto_F$adj.P.Val <0.05] <- "Trending Towards Increased Expression"
affy_Lip_auto_F$Expression_Status[affy_Lip_auto_F$logFC < 0.00 & affy_Lip_auto_F$adj.P.Val <0.05] <- "Trending Towards Decreased Expression"
affy_Lip_auto_F$Expression_Status[affy_Lip_auto_F$logFC > 0.05 & affy_Lip_auto_F$adj.P.Val <0.05] <- "Increased Expression"
affy_Lip_auto_F$Expression_Status[affy_Lip_auto_F$logFC < -0.05 & affy_Lip_auto_F$adj.P.Val <0.05] <- "Decreased Expression"

#Create List of significant and biologically significant genes for each platform

# Female 
DNAm_sig_F <- unique(DNAm_Lip_auto_F$gene[DNAm_Lip_auto_F$diffmethylation != "Not_Biologically_Significant"])
DNAm_biosig_F <- unique(DNAm_Lip_auto_F$gene[DNAm_Lip_auto_F$diffmethylation %in% c("Increased Methylation", "Decreased Methylation")])
RNA_sig_F <- unique(RNA_Lip_auto_F$gene_symbol[RNA_Lip_auto_F$Expression_Status != "Not_Biologically_Significant"])
RNA_biosig_F <- unique(RNA_Lip_auto_F$gene_symbol[RNA_Lip_auto_F$Expression_Status %in% c("Increased_RNA_Expression", "Decreased_RNA_Expression")])
affy_sig_F <- unique(affy_Lip_auto_F$SYMBOL[affy_Lip_auto_F$Expression_Status != "Not_Biologically_Significant"])
affy_biosig_F <- unique(affy_Lip_auto_F$SYMBOL[affy_Lip_auto_F$Expression_Status %in% c("Increased Expression", "Decreased Expression")])

#Create list of genes for each category in upset plot
female_upset_list <- list("DNAm_sig_F" = DNAm_sig_F,"DNAm_biosig_F" = DNAm_biosig_F, "RNA_sig_F" = RNA_sig_F, "RNA_biosig_F" = RNA_biosig_F, "Affy_sig_F" = affy_sig_F, "Affy_biosig_F" = affy_biosig_F)

order_F <- c("DNAm_biosig_F", "DNAm_sig_F", "Affy_biosig_F", "Affy_sig_F", "RNA_biosig_F", "RNA_sig_F")
png(file = "./Female_Tiered_Overlap_UpSet.png", height = 8, width = 12, units = "in", res = 300)
upset(fromList(female_upset_list), 
      sets = order_F,
      keep.order = TRUE,
      mainbar.y.label = "Number of Shared Genes",
      nsets = 6, 
      order.by = "freq", 
      sets.bar.color = "black", 
      main.bar.color = "#d02670",
      text.scale = 1.75)
dev.off()

#Identify which gene is overlapping between the different platforms 
f_upset_table <- as.data.frame(fromList(female_upset_list))
all_gene_names <- unique(unlist(female_upset_list))
rownames(f_upset_table) <- all_gene_names
head(f_upset_table) #RAB5C

f_upset_table$overlap_count <- rowSums(f_upset_table)
table(f_upset_table$overlap_count)
f_upset_table[f_upset_table$overlap_count == 4, ] #NA



# Male 
DNAm_sig_M <- unique(DNAm_Lip_auto_M$gene[DNAm_Lip_auto_M$diffmethylation != "Not_Biologically_Significant"])
DNAm_biosig_M <- unique(DNAm_Lip_auto_M$gene[DNAm_Lip_auto_M$diffmethylation %in% c("Increased Methylation", "Decreased Methylation")])
RNA_sig_M <- unique(RNA_Lip_auto_M$gene_symbol[RNA_Lip_auto_M$Expression_Status != "Not_Biologically_Significant"])
RNA_biosig_M <- unique(RNA_Lip_auto_M$gene_symbol[RNA_Lip_auto_M$Expression_Status %in% c("Increased_RNA_Expression", "Decreased_RNA_Expression")])
affy_sig_M <- unique(affy_Lip_auto_M$SYMBOL[affy_Lip_auto_M$Expression_Status != "Not_Biologically_Significant"])
affy_biosig_M <- unique(affy_Lip_auto_M$SYMBOL[affy_Lip_auto_M$Expression_Status %in% c("Increased Expression", "Decreased Expression")])

DNAm_sig_M <- DNAm_sig_M[!is.na(DNAm_sig_M)]
DNAm_biosig_M <- DNAm_biosig_M[!is.na(DNAm_biosig_M)]
RNA_sig_M <- RNA_sig_M[!is.na(RNA_sig_M)]
RNA_biosig_M <- RNA_biosig_M[!is.na(RNA_biosig_M)]
affy_sig_M <- affy_sig_M[!is.na(affy_sig_M)]
affy_biosig_M <- affy_biosig_M[!is.na(affy_biosig_M)]

#Create list of genes for each category in upset plot

male_upset_list <- list("DNAm_sig_M" = DNAm_sig_M,"DNAm_biosig_M" = DNAm_biosig_M, "RNA_sig_M" = RNA_sig_M, "RNA_biosig_M" = RNA_biosig_M, "Affy_sig_M" = affy_sig_M, "Affy_biosig_M" = affy_biosig_M)
male_upset_list <- male_upset_list[!is.na(male_upset_list)]

order_M <- c("DNAm_biosig_M", "DNAm_sig_M", "Affy_biosig_M", "Affy_sig_M", "RNA_biosig_M", "RNA_sig_M")
png(file = "./Male_Tiered_Overlap_UpSet.png", height = 8, width = 12, units = "in", res = 300)
upset(fromList(male_upset_list), 
      sets = order_M,
      keep.order = TRUE,
      mainbar.y.label = "Number of Shared Genes",
      nsets = 6, 
      order.by = "freq", 
      sets.bar.color = "black", 
      main.bar.color = "#275317",
      text.scale = 1.75)
dev.off()

#Identify which gene is overlapping between the different platforms 
m_upset_table <- as.data.frame(fromList(male_upset_list))
all_gene_names_m <- unique(unlist(male_upset_list))
rownames(m_upset_table) <- all_gene_names_m
head(m_upset_table) #GATA3 - RNABiosig, SH3PXD2A, SLC2A1

m_upset_table$overlap_count <- rowSums(m_upset_table)
table(m_upset_table$overlap_count)
m_upset_table[m_upset_table$overlap_count == 5, ] #SH3PXD2A, C8orf58


#Dual Overlap Upset Plot
# put image in ggplot2 object
library(grid)
library(gridExtra)
library(png)
img_f <- rasterGrob(readPNG("Female_Tiered_Overlap_UpSet.png"), interpolate=TRUE)
img_m <- rasterGrob(readPNG("Male_Tiered_Overlap_UpSet.png"), interpolate=TRUE)

png("./Female_Male_Overlap_Upset.png", height = 9, width = 20, units = "in", res = 300)
grid.arrange(img_f, img_m, nrow = 1)
dev.off()



# Upset Plot looking at overlap of differentially expressed genes (RNA-seq) between males and females
RNA_sig_increased_F <- unique(RNA_Lip_auto_F$gene_symbol[RNA_Lip_auto_F$Expression_Status %in% c("Increased_RNA_Expression", "Trending_Towards_Increased_RNA_Expression")])
RNA_sig_decreased_F <- unique(RNA_Lip_auto_F$gene_symbol[RNA_Lip_auto_F$Expression_Status %in% c("Decreased_RNA_Expression", "Trending_Towards_Decreased_RNA_Expression")])
RNA_biosig_increased_F <- unique(RNA_Lip_auto_F$gene_symbol[RNA_Lip_auto_F$Expression_Status %in% c("Increased_RNA_Expression")])
RNA_biosig_decreased_F <- unique(RNA_Lip_auto_F$gene_symbol[RNA_Lip_auto_F$Expression_Status %in% c("Decreased_RNA_Expression")])

RNA_sig_increased_M <- unique(RNA_Lip_auto_M$gene_symbol[RNA_Lip_auto_M$Expression_Status %in% c("Increased_RNA_Expression", "Trending_Towards_Increased_RNA_Expression")])
RNA_sig_decreased_M <- unique(RNA_Lip_auto_M$gene_symbol[RNA_Lip_auto_M$Expression_Status %in% c("Decreased_RNA_Expression", "Trending_Towards_Decreased_RNA_Expression")])
RNA_biosig_increased_M <- unique(RNA_Lip_auto_M$gene_symbol[RNA_Lip_auto_M$Expression_Status %in% c("Increased_RNA_Expression")])
RNA_biosig_decreased_M <- unique(RNA_Lip_auto_M$gene_symbol[RNA_Lip_auto_M$Expression_Status %in% c("Decreased_RNA_Expression")])

female_male_overlap_list <- list("RNA_sig_increased_F" = RNA_sig_increased_F,"RNA_sig_decreased_F" = RNA_sig_decreased_F, "RNA_biosig_increased_F" = RNA_biosig_increased_F, "RNA_biosig_decreased_F" = RNA_biosig_decreased_F, "RNA_sig_increased_M" = RNA_sig_increased_M, "RNA_sig_decreased_M" = RNA_sig_decreased_M, "RNA_biosig_increased_M" = RNA_biosig_increased_M, "RNA_biosig_decreased_M" = RNA_biosig_decreased_M)
female_male_biosig_overlap_list <- list("RNA_biosig_increased_F" = RNA_biosig_increased_F, "RNA_biosig_decreased_F" = RNA_biosig_decreased_F, "RNA_biosig_increased_M" = RNA_biosig_increased_M, "RNA_biosig_decreased_M" = RNA_biosig_decreased_M)


order <- c("RNA_sig_increased_F", "RNA_biosig_increased_F", "RNA_sig_decreased_F", "RNA_biosig_decreased_F", "RNA_sig_increased_M", "RNA_biosig_increased_M", "RNA_sig_decreased_M", "RNA_biosig_decreased_M")
png(file = "./Female_Male_RNA_Overlap_UpSet.png", height = 8, width = 12, units = "in", res = 300)
upset(fromList(female_male_overlap_list), 
      sets = order,
      keep.order = TRUE,
      mainbar.y.label = "Number of Shared Genes",
      nsets = 8, 
      order.by = "freq", 
      sets.bar.color = "black", 
      main.bar.color = "black",
      text.scale = 1.75)
dev.off()

order_biosig <- c("RNA_biosig_increased_F", "RNA_biosig_decreased_F", "RNA_biosig_increased_M", "RNA_biosig_decreased_M")
png(file = "./Female_Male_RNA_Biosig_Overlap_UpSet.png", height = 8, width = 12, units = "in", res = 300)
upset(fromList(female_male_biosig_overlap_list), 
      sets = order_biosig,
      keep.order = TRUE,
      mainbar.y.label = "Number of Shared Genes",
      nsets = 8, 
      order.by = "freq", 
      sets.bar.color = "black", 
      main.bar.color = "black",
      text.scale = 1.75)
dev.off()
biosig_upset_table <- as.data.frame(fromList(female_male_biosig_overlap_list))
all_gene_names_biosig <- unique(unlist(female_male_biosig_overlap_list))
rownames(biosig_upset_table) <- all_gene_names_biosig
head(biosig_upset_table) #NR4A1

decreased_M <- biosig_upset_table[biosig_upset_table$RNA_biosig_decreased_M =="1",]
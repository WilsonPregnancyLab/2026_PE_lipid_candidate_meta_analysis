# Packages (Run)
library(dplyr) #version 1.1.4


##Load Candidate Differential RNA Expression Results
RNA_Lip_DESeq_results_autosomes_combined_sex <- read.csv("/workspace/lab/wilsonslab/datalake-wilsonslab/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/DESeq_analysis/RNA_Lip_DESeq_results_autosomes_combined_sex.csv")
RNA_Lip_DESeq_results_autosomes_F <- read.csv("/workspace/lab/wilsonslab/datalake-wilsonslab/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/DESeq_analysis/RNA_Lip_DESeq_results_autosomes_F.csv")
RNA_Lip_DESeq_results_autosomes_M <- read.csv("/workspace/lab/wilsonslab/datalake-wilsonslab/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/DESeq_analysis/RNA_Lip_DESeq_results_autosomes_M.csv")
RNA_Lip_DESeq_results_chrX_F <- read.csv("/workspace/lab/wilsonslab/datalake-wilsonslab/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/DESeq_analysis/RNA_Lip_DESeq_results_chrX_F.csv")
RNA_Lip_DESeq_results_chrX_M <- read.csv("/workspace/lab/wilsonslab/datalake-wilsonslab/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/DESeq_analysis/RNA_Lip_DESeq_results_chrX_M.csv")


##Load Candidate Differential DNA Methylation Results

lipid_candidate_probes <- read.csv("/workspace/lab/wilsonslab/datalake-wilsonslab/PE_Lipid_Meta-analysis/2025_DNAm_Lipid_Candidate/Candidate_Diff_DNAm_Analysis_2025/lipid_candidate_probes.csv")
placmet_wholepop_auto <- read.csv ("/workspace/lab/wilsonslab/datalake-wilsonslab/PE_Lipid_Meta-analysis/2025_DNAm_Lipid_Candidate/Candidate_Diff_DNAm_Analysis_2025/placmet_wholepop_auto_rerun.csv")
placmet_M_fulldata_auto <- read.csv("/workspace/lab/wilsonslab/datalake-wilsonslab/PE_Lipid_Meta-analysis/2025_DNAm_Lipid_Candidate/Candidate_Diff_DNAm_Analysis_2025/placmet_M_fulldata_auto_rerun.csv")
placmet_M_fulldata_X <- read.csv("/workspace/lab/wilsonslab/datalake-wilsonslab/PE_Lipid_Meta-analysis/2025_DNAm_Lipid_Candidate/Candidate_Diff_DNAm_Analysis_2025/placmet_M_fulldata_X_rerun.csv")
placmet_F_fulldata_auto <- read.csv("/workspace/lab/wilsonslab/datalake-wilsonslab/PE_Lipid_Meta-analysis/2025_DNAm_Lipid_Candidate/Candidate_Diff_DNAm_Analysis_2025/placmet_F_fulldata_auto_rerun.csv")
placmet_F_fulldata_X <- read.csv("/workspace/lab/wilsonslab/datalake-wilsonslab/PE_Lipid_Meta-analysis/2025_DNAm_Lipid_Candidate/Candidate_Diff_DNAm_Analysis_2025/placmet_F_fulldata_X_rerun.csv")

plot_wholepop_auto <- placmet_wholepop_auto[placmet_wholepop_auto$gene %in% lipid_candidate_probes$gene,]
plot_M_fulldata_auto <- placmet_M_fulldata_auto[placmet_M_fulldata_auto$gene %in% lipid_candidate_probes$gene,]
plot_M_fulldata_X <- placmet_M_fulldata_X[placmet_M_fulldata_X$gene %in% lipid_candidate_probes$gene,]
plot_F_fulldata_auto <- placmet_F_fulldata_auto[placmet_F_fulldata_auto$gene %in% lipid_candidate_probes$gene,]
plot_F_fulldata_X <- placmet_F_fulldata_X[placmet_F_fulldata_X$gene %in% lipid_candidate_probes$gene,]

#Whole data 
plot_wholepop_auto$diffmethylation <- "Not_Biologically_Significant"
plot_wholepop_auto$diffmethylation[plot_wholepop_auto$deltaB > 0.00 & plot_wholepop_auto$adj.P.Val <0.05] <- "Trending Towards Increased Methylation"
plot_wholepop_auto$diffmethylation[plot_wholepop_auto$deltaB < 0.00 & plot_wholepop_auto$adj.P.Val <0.05] <- "Trending Towards Decreased Methylation"
plot_wholepop_auto$diffmethylation[plot_wholepop_auto$deltaB > 0.05 & plot_wholepop_auto$adj.P.Val <0.05] <- "Increased Methylation"
plot_wholepop_auto$diffmethylation[plot_wholepop_auto$deltaB < -0.05 & plot_wholepop_auto$adj.P.Val <0.05] <- "Decreased Methylation"
#Male data 
plot_M_fulldata_auto$diffmethylation <- "Not_Biologically_Significant"
plot_M_fulldata_auto$diffmethylation[plot_M_fulldata_auto$deltaB > 0.00 & plot_M_fulldata_auto$adj.P.Val <0.05] <- "Trending Towards Increased Methylation"
plot_M_fulldata_auto$diffmethylation[plot_M_fulldata_auto$deltaB < 0.00 & plot_M_fulldata_auto$adj.P.Val <0.05] <- "Trending Towards Decreased Methylation"
plot_M_fulldata_auto$diffmethylation[plot_M_fulldata_auto$deltaB > 0.05 & plot_M_fulldata_auto$adj.P.Val <0.05] <- "Increased Methylation"
plot_M_fulldata_auto$diffmethylation[plot_M_fulldata_auto$deltaB < -0.05 & plot_M_fulldata_auto$adj.P.Val <0.05] <- "Decreased Methylation"
plot_M_fulldata_X$diffmethylation <- "Not_Biologically_Significant"
plot_M_fulldata_X$diffmethylation[plot_M_fulldata_X$deltaB > 0.00 & plot_M_fulldata_X$adj.P.Val <0.05] <- "Trending Towards Increased Methylation"
plot_M_fulldata_X$diffmethylation[plot_M_fulldata_X$deltaB < 0.00 & plot_M_fulldata_X$adj.P.Val <0.05] <- "Trending Towards Decreased Methylation"
#Female Data 
plot_F_fulldata_auto$diffmethylation <- "Not_Biologically_Significant"
plot_F_fulldata_auto$diffmethylation[plot_F_fulldata_auto$deltaB > 0.00 & plot_F_fulldata_auto$adj.P.Val <0.05] <- "Trending Towards Increased Methylation"
plot_F_fulldata_auto$diffmethylation[plot_F_fulldata_auto$deltaB < 0.00 & plot_F_fulldata_auto$adj.P.Val <0.05] <- "Trending Towards Decreased Methylation"
plot_F_fulldata_auto$diffmethylation[plot_F_fulldata_auto$deltaB > 0.05 & plot_F_fulldata_auto$adj.P.Val <0.05] <- "Increased Methylation"
plot_F_fulldata_auto$diffmethylation[plot_F_fulldata_auto$deltaB < -0.05 & plot_F_fulldata_auto$adj.P.Val <0.05] <- "Decreased Methylation"
plot_F_fulldata_X$diffmethylation <- "Not_Biologically_Significant"
plot_F_fulldata_X$diffmethylation[plot_F_fulldata_X$deltaB > 0.00 & plot_F_fulldata_X$adj.P.Val <0.05] <- "Trending Towards Increased Methylation"
plot_F_fulldata_X$diffmethylation[plot_F_fulldata_X$deltaB < 0.00 & plot_F_fulldata_X$adj.P.Val <0.05] <- "Trending Towards Decreased Methylation"


wholeauto_met_sig <- subset(plot_wholepop_auto[plot_wholepop_auto$adj.P.Val <0.05,])
femaleauto_met_sig <- subset(plot_F_fulldata_auto[plot_F_fulldata_auto$adj.P.Val <0.05,]) #2
maleauto_met_sig <- subset(plot_M_fulldata_auto[plot_M_fulldata_auto$adj.P.Val <0.05,]) #16
wholeauto_met_bio_sig<- subset(plot_wholepop_auto[plot_wholepop_auto$adj.P.Val <0.05 & (plot_wholepop_auto$deltaB < -0.05 | plot_wholepop_auto$deltaB > 0.05), ]) #0
femaleauto_met_bio_sig <- subset(plot_F_fulldata_auto[plot_F_fulldata_auto$adj.P.Val <0.05 & (plot_F_fulldata_auto$deltaB < -0.05 | plot_F_fulldata_auto$deltaB > 0.05), ]) #2
maleauto_met_bio_sig <- subset(plot_M_fulldata_auto[plot_M_fulldata_auto$adj.P.Val <0.05 & (plot_M_fulldata_auto$deltaB < -0.05 | plot_M_fulldata_auto$deltaB > 0.05), ]) #10

#Joined Table

##Autosomes Female
F_RNA_auto_subset <- RNA_Lip_DESeq_results_autosomes_F[,c("Row.names", "log2FoldChange", "padj", "gene_symbol", "seqnames", "Expression_Status")]
F_met_auto_subset <- femaleauto_met_sig[,c("probes", "adj.P.Val", "deltaB", "gene", "chr", "position", "region_overlap", "Closest_TSS_gene_name", "Closest_TSS_Transcript", "Closest_TSS_Pos", "Distance_Closest_TSS", "gene_id", "overlap", "diffmethylation")]

overlap_F <- F_RNA_auto_subset[F_RNA_auto_subset$gene_symbol %in% femaleauto_met_sig$gene, ]
auto_F_table <- right_join(F_RNA_auto_subset, F_met_auto_subset, by = c("gene_symbol" = "gene"))
write.csv(auto_F_table, file = "Lip_PE_Met_RNA_auto_F.csv")

##Autosomes Male
M_RNA_auto_subset <- RNA_Lip_DESeq_results_autosomes_M[,c("Row.names", "log2FoldChange", "padj", "gene_symbol", "seqnames", "Expression_Status")]
M_met_auto_subset <- maleauto_met_sig[,c("probes", "adj.P.Val", "deltaB", "gene", "chr", "position", "region_overlap", "Closest_TSS_gene_name", "Closest_TSS_Transcript", "Closest_TSS_Pos", "Distance_Closest_TSS", "gene_id", "overlap", "diffmethylation")]

overlap_M <- M_RNA_auto_subset[M_RNA_auto_subset$gene_symbol %in% maleauto_met_sig$gene, ]
auto_M_table <- right_join(M_RNA_auto_subset, M_met_auto_subset, by = c("gene_symbol" = "gene"))
write.csv(auto_M_table, file = "Lip_PE_Met_RNA_auto_M.csv")

Lip_PE_Met_RNA_auto_table <- full_join(auto_F_table, auto_M_table, by = c("Row.names", "gene_symbol", "seqnames", "probes", "chr", "position", "region_overlap", "Closest_TSS_gene_name", "Closest_TSS_Transcript", "Closest_TSS_Pos", "Distance_Closest_TSS", "gene_id", "overlap"))
write.csv(Lip_PE_Met_RNA_auto_table, "Lip_PE_Met_RNA_auto_table.csv")


#overlap in affymetrix validation cohort
affy_results_auto_F <- read.csv("/workspace/lab/wilsonslab/datalake-wilsonslab/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/affymetrix_validation/affy_results_auto_F.csv")
affy_results_auto_M <- read.csv("/workspace/lab/wilsonslab/datalake-wilsonslab/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/affymetrix_validation/affy_results_auto_M.csv")
female_auto_sig <- subset(affy_results_auto_F[affy_results_auto_F$adj.P.Val <0.05,]) # 0 
male_auto_sig <- subset(affy_results_auto_M[affy_results_auto_M$adj.P.Val <0.05,]) # 0 

affy_overlap_F <- affy_results_auto_F[affy_results_auto_F$SYMBOL %in% female_auto_sig$gene,] #no sig differentially methylated CpGs or genes overlap with affymetrix cohort
affy_overlap_M <- affy_results_auto_M[affy_results_auto_M$SYMBOL %in% male_auto_sig$gene,] #no sig differentially methylated CpGs or genes overlap with affymetrix cohort






##WG Overlap

##Load WG Differential RNA Expression Results
RNA_WG_DESeq_results_autosomes_combined_sex <- read.csv("/workspace/lab/wilsonslab/datalake-wilsonslab/PE_Lipid_Meta-analysis/2025_RNA_WG/RNA_WG_DESeq_results_autosomes_combined_sex.csv")
RNA_WG_DESeq_results_autosomes_F <- read.csv("/workspace/lab/wilsonslab/datalake-wilsonslab/PE_Lipid_Meta-analysis/2025_RNA_WG/RNA_WG_DESeq_results_autosomes_F.csv")
RNA_WG_DESeq_results_autosomes_M <- read.csv("/workspace/lab/wilsonslab/datalake-wilsonslab/PE_Lipid_Meta-analysis/2025_RNA_WG/RNA_WG_DESeq_results_autosomes_M.csv")
RNA_WG_DESeq_results_chrX_F <- read.csv("/workspace/lab/wilsonslab/datalake-wilsonslab/PE_Lipid_Meta-analysis/2025_RNA_WG/RNA_WG_DESeq_results_chrX_F.csv")
RNA_WG_DESeq_results_chrX_M <- read.csv("/workspace/lab/wilsonslab/datalake-wilsonslab/PE_Lipid_Meta-analysis/2025_RNA_WG/RNA_WG_DESeq_results_chrX_M.csv")
RNA_WG_DESeq_results_chrY_M <- read.csv("/workspace/lab/wilsonslab/datalake-wilsonslab/PE_Lipid_Meta-analysis/2025_RNA_WG/RNA_WG_DESeq_results_chrY_M.csv")

##Load WG Differential DNA Methylation Results

WG_placmet_wholepop_auto <- read.csv ("/workspace/lab/wilsonslab/datalake-wilsonslab/PE_Lipid_Meta-analysis/2025_DNAm_Lipid_Candidate/WG_Diff_DNAm_Analysis_2025/wg_placmet_wholepop_auto_rerun.csv")
WG_placmet_F_fulldata_auto <- read.csv("/workspace/lab/wilsonslab/datalake-wilsonslab/PE_Lipid_Meta-analysis/2025_DNAm_Lipid_Candidate/WG_Diff_DNAm_Analysis_2025/wg_placmet_F_fulldata_auto_rerun.csv")
WG_placmet_F_fulldata_X <- read.csv("/workspace/lab/wilsonslab/datalake-wilsonslab/PE_Lipid_Meta-analysis/2025_DNAm_Lipid_Candidate/WG_Diff_DNAm_Analysis_2025/wg_placmet_F_fulldata_X_rerun.csv")
WG_placmet_M_fulldata_auto <- read.csv("/workspace/lab/wilsonslab/datalake-wilsonslab/PE_Lipid_Meta-analysis/2025_DNAm_Lipid_Candidate/WG_Diff_DNAm_Analysis_2025/wg_placmet_M_fulldata_auto_rerun.csv")
WG_placmet_M_fulldata_X <- read.csv("/workspace/lab/wilsonslab/datalake-wilsonslab/PE_Lipid_Meta-analysis/2025_DNAm_Lipid_Candidate/WG_Diff_DNAm_Analysis_2025/wg_placmet_M_fulldata_X_rerun.csv")
WG_placmet_M_fulldata_Y <- read.csv("/workspace/lab/wilsonslab/datalake-wilsonslab/PE_Lipid_Meta-analysis/2025_DNAm_Lipid_Candidate/WG_Diff_DNAm_Analysis_2025/wg_placmet_M_fulldata_Y_rerun.csv")

#Whole data 
WG_placmet_wholepop_auto$diffmethylation <- "Not_Biologically_Significant"
WG_placmet_wholepop_auto$diffmethylation[WG_placmet_wholepop_auto$deltaB > 0.00 & WG_placmet_wholepop_auto$adj.P.Val <0.05] <- "Trending Towards Increased Methylation"
WG_placmet_wholepop_auto$diffmethylation[WG_placmet_wholepop_auto$deltaB < 0.00 & WG_placmet_wholepop_auto$adj.P.Val <0.05] <- "Trending Towards Decreased Methylation"
WG_placmet_wholepop_auto$diffmethylation[WG_placmet_wholepop_auto$deltaB > 0.05 & WG_placmet_wholepop_auto$adj.P.Val <0.05] <- "Increased Methylation"
WG_placmet_wholepop_auto$diffmethylation[WG_placmet_wholepop_auto$deltaB < -0.05 & WG_placmet_wholepop_auto$adj.P.Val <0.05] <- "Decreased Methylation"
#Male data 
WG_placmet_M_fulldata_auto$diffmethylation <- "Not_Biologically_Significant"
WG_placmet_M_fulldata_auto$diffmethylation[WG_placmet_M_fulldata_auto$deltaB > 0.00 & WG_placmet_M_fulldata_auto$adj.P.Val <0.05] <- "Trending Towards Increased Methylation"
WG_placmet_M_fulldata_auto$diffmethylation[WG_placmet_M_fulldata_auto$deltaB < 0.00 & WG_placmet_M_fulldata_auto$adj.P.Val <0.05] <- "Trending Towards Decreased Methylation"
WG_placmet_M_fulldata_auto$diffmethylation[WG_placmet_M_fulldata_auto$deltaB > 0.05 & WG_placmet_M_fulldata_auto$adj.P.Val <0.05] <- "Increased Methylation"
WG_placmet_M_fulldata_auto$diffmethylation[WG_placmet_M_fulldata_auto$deltaB < -0.05 & WG_placmet_M_fulldata_auto$adj.P.Val <0.05] <- "Decreased Methylation"
WG_placmet_M_fulldata_X$diffmethylation <- "Not_Biologically_Significant"
WG_placmet_M_fulldata_X$diffmethylation[WG_placmet_M_fulldata_X$deltaB > 0.00 & WG_placmet_M_fulldata_X$adj.P.Val <0.05] <- "Trending Towards Increased Methylation"
WG_placmet_M_fulldata_X$diffmethylation[WG_placmet_M_fulldata_X$deltaB < 0.00 & WG_placmet_M_fulldata_X$adj.P.Val <0.05] <- "Trending Towards Decreased Methylation"
#Female Data 
WG_placmet_F_fulldata_auto$diffmethylation <- "Not_Biologically_Significant"
WG_placmet_F_fulldata_auto$diffmethylation[WG_placmet_F_fulldata_auto$deltaB > 0.00 & WG_placmet_F_fulldata_auto$adj.P.Val <0.05] <- "Trending Towards Increased Methylation"
WG_placmet_F_fulldata_auto$diffmethylation[WG_placmet_F_fulldata_auto$deltaB < 0.00 & WG_placmet_F_fulldata_auto$adj.P.Val <0.05] <- "Trending Towards Decreased Methylation"
WG_placmet_F_fulldata_auto$diffmethylation[WG_placmet_F_fulldata_auto$deltaB > 0.05 & WG_placmet_F_fulldata_auto$adj.P.Val <0.05] <- "Increased Methylation"
WG_placmet_F_fulldata_auto$diffmethylation[WG_placmet_F_fulldata_auto$deltaB < -0.05 & WG_placmet_F_fulldata_auto$adj.P.Val <0.05] <- "Decreased Methylation"

WG_placmet_F_fulldata_X$diffmethylation <- "Not_Biologically_Significant"
WG_placmet_F_fulldata_X$diffmethylation[WG_placmet_F_fulldata_X$deltaB > 0.00 & WG_placmet_F_fulldata_X$adj.P.Val <0.05] <- "Trending Towards Increased Methylation"
WG_placmet_F_fulldata_X$diffmethylation[WG_placmet_F_fulldata_X$deltaB < 0.00 & WG_placmet_F_fulldata_X$adj.P.Val <0.05] <- "Trending Towards Decreased Methylation"

WG_wholeauto_met_sig <- subset(WG_placmet_wholepop_auto[WG_placmet_wholepop_auto$adj.P.Val <0.05,]) #0
WG_femaleauto_met_sig <- subset(WG_placmet_F_fulldata_auto[WG_placmet_F_fulldata_auto$adj.P.Val <0.05,]) #2
WG_maleauto_met_sig <- subset(WG_placmet_M_fulldata_auto[WG_placmet_M_fulldata_auto$adj.P.Val <0.05,]) #16
WG_femalechrX_met_sig <- subset(WG_placmet_F_fulldata_X[WG_placmet_F_fulldata_X$adj.P.Val <0.05,]) #0
WG_malechrX_met_sig <- subset(WG_placmet_M_fulldata_X[WG_placmet_M_fulldata_X$adj.P.Val <0.05,]) #0
WG_malechrY_met_sig <- subset(WG_placmet_M_fulldata_Y[WG_placmet_M_fulldata_Y$adj.P.Val <0.05,]) #0

WG_femaleauto_met_bio_sig <- subset(WG_placmet_F_fulldata_auto[WG_placmet_F_fulldata_auto$adj.P.Val <0.05 & (WG_placmet_F_fulldata_auto$deltaB < -0.05 | WG_placmet_F_fulldata_auto$deltaB > 0.05), ]) #2
WG_maleauto_met_bio_sig <- subset(WG_placmet_M_fulldata_auto[WG_placmet_M_fulldata_auto$adj.P.Val <0.05 & (WG_placmet_M_fulldata_auto$deltaB < -0.05 | WG_placmet_M_fulldata_auto$deltaB > 0.05), ]) #10

#Joined Table

##Autosomes Female
F_RNA_auto_subset_WG <- RNA_WG_DESeq_results_autosomes_F[,c("Row.names", "log2FoldChange", "padj", "gene_name", "seqnames", "Expression_Status")]
F_met_auto_subset_WG <- WG_femaleauto_met_sig[,c("probes", "adj.P.Val", "deltaB", "gene", "chr", "position", "region_overlap", "Closest_TSS_gene_name", "Closest_TSS_Transcript", "Closest_TSS_Pos", "Distance_Closest_TSS", "gene_id", "overlap", "diffmethylation")]

WG_overlap_F <- F_RNA_auto_subset_WG[F_RNA_auto_subset_WG$gene_name %in% WG_femaleauto_met_sig$gene, ] #0
WG_auto_F_table <- right_join(F_RNA_auto_subset_WG, F_met_auto_subset_WG, by = c("gene_name" = "gene"))
write.csv(WG_auto_F_table, file = "WG_PE_Met_RNA_auto_F.csv")

##Autosomes Male
M_RNA_auto_subset_WG <- RNA_WG_DESeq_results_autosomes_M[,c("Row.names", "log2FoldChange", "padj", "gene_name", "seqnames", "Expression_Status")]
M_met_auto_subset_WG <- WG_maleauto_met_sig[,c("probes", "adj.P.Val", "deltaB", "gene", "chr", "position", "region_overlap", "Closest_TSS_gene_name", "Closest_TSS_Transcript", "Closest_TSS_Pos", "Distance_Closest_TSS", "gene_id", "overlap", "diffmethylation")]

WG_overlap_M <- M_RNA_auto_subset_WG[M_RNA_auto_subset_WG$gene_name %in% WG_maleauto_met_sig$gene, ]
WG_auto_M_table <- right_join(M_RNA_auto_subset_WG, M_met_auto_subset_WG, by = c("gene_name" = "gene"))
write.csv(WG_auto_M_table, file = "WG_PE_Met_RNA_auto_M.csv")

WG_PE_Met_RNA_auto_table <- full_join(WG_auto_F_table, WG_auto_M_table, by = c("Row.names", "gene_name", "seqnames", "probes", "chr", "position", "region_overlap", "Closest_TSS_gene_name", "Closest_TSS_Transcript", "Closest_TSS_Pos", "Distance_Closest_TSS", "gene_id", "overlap"))
write.csv(WG_PE_Met_RNA_auto_table, "WG_PE_Met_RNA_auto_table.csv")

#overlap in affymetrix validation cohort
WG_affy_results_auto_F <- read.csv("/workspace/lab/wilsonslab/datalake-wilsonslab/PE_Lipid_Meta-analysis/2025_RNA_WG/affymetrix_validation_WG/affy_results_auto_F.csv")
WG_affy_results_auto_M <- read.csv("/workspace/lab/wilsonslab/datalake-wilsonslab/PE_Lipid_Meta-analysis/2025_RNA_WG/affymetrix_validation_WG/affy_results_auto_M.csv")
female_auto_sig <- subset(WG_affy_results_auto_F[WG_affy_results_auto_F$adj.P.Val <0.05,]) # 0 
male_auto_sig <- subset(WG_affy_results_auto_M[WG_affy_results_auto_M$adj.P.Val <0.05,]) # 0 

WG_affy_overlap_F <- WG_affy_results_auto_F[WG_affy_results_auto_F$SYMBOL %in% female_auto_sig$gene,] #no sig differentially methylated CpGs or genes overlap with affymetrix cohort
WG_affy_overlap_M <- WG_affy_results_auto_M[WG_affy_results_auto_M$SYMBOL %in% male_auto_sig$gene,] #no sig differentially methylated CpGs or genes overlap with affymetrix cohort





















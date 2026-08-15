# Comparison Table
library(dplyr) #version 1.1.4


Overlap_Chart_Lip <- function(DNAm_file, RNA_seq_file, affy_file){
    affy_name <- deparse(substitute(affy_file))
    no_affy_name <- gsub("affy_", "", affy_name)
    filename <- paste0("overlap_", no_affy_name, ".csv")
    print(filename)
    
    DNAm <- DNAm_file[, colnames(DNAm_file) %in% c("gene", "Closest_TSS_gene_name", "probes", "gene_id", "logFC", "adj.P.Val", "deltaB", "region_overlap", "chr")]
    sig_DNAm <- DNAm[DNAm$adj.P.Val < 0.05,]
    RNA_Seq <- RNA_seq_file[, colnames(RNA_seq_file) %in% c("Row.names","gene_symbol", "log2FoldChange", "padj", "Expression_Status")]
    sig_RNA_Seq <- RNA_Seq[RNA_Seq$padj < 0.05,]
    affy <- affy_file[, colnames(affy_file) %in% c("ENSEMBL", "SYMBOL", "adj.P.Val", "logFC", "deltaExprs")]

    RNAmerge <- merge(sig_RNA_Seq, affy, by.x = "gene_symbol", by.y = "SYMBOL", all = TRUE)
    overlap_output <- merge(RNAmerge, DNAm, by.x = "gene_symbol", by.y = "gene", all = TRUE)
    overlap_output <- merge(overlap_output, gtf_df[ , c("gene_id", "gene_type")], by = "gene_id")

    write.csv(overlap_output, file = filename)
}

DNAm_Lip_auto_combinedsex <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_DNAm_Lipid_Candidate/Candidate_Diff_DNAm_Analysis_2025/placmet_wholepop_auto_rerun.csv")
DNAm_Lip_auto_F <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_DNAm_Lipid_Candidate/Candidate_Diff_DNAm_Analysis_2025/placmet_F_fulldata_auto_rerun.csv")
DNAm_Lip_auto_M <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_DNAm_Lipid_Candidate/Candidate_Diff_DNAm_Analysis_2025/placmet_M_fulldata_auto_rerun.csv")
DNAm_Lip_chrX_F <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_DNAm_Lipid_Candidate/Candidate_Diff_DNAm_Analysis_2025/placmet_F_fulldata_X_rerun.csv")
DNAm_Lip_chrX_M <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_DNAm_Lipid_Candidate/Candidate_Diff_DNAm_Analysis_2025/placmet_M_fulldata_X_rerun.csv")

RNA_Lip_auto_combinedsex <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/DESeq_rerun/RNA_Lip_DESeq_results_autosomes_combined_sex.csv")
RNA_Lip_auto_F <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/DESeq_rerun/RNA_Lip_DESeq_results_autosomes_F.csv")
RNA_Lip_auto_M <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/DESeq_rerun/RNA_Lip_DESeq_results_autosomes_M.csv")
RNA_Lip_chrX_F <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/DESeq_rerun/RNA_Lip_DESeq_results_chrX_F.csv")
RNA_Lip_chrX_M <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/DESeq_rerun/RNA_Lip_DESeq_results_chrX_M.csv")

affy_Lip_auto_combinedsex <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/affymetrix_validation_rerun/affy_results_auto_combined_sex.csv") 
affy_Lip_auto_F <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/affymetrix_validation_rerun/affy_results_auto_F.csv")
affy_Lip_auto_M <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/affymetrix_validation_rerun/affy_results_auto_M.csv")
affy_Lip_chrX_F <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/affymetrix_validation_rerun/affy_results_chrX_F.csv")
affy_Lip_chrX_M <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/affymetrix_validation_rerun/affy_results_chrX_M.csv")

Overlap_Chart_Lip(DNAm_Lip_auto_combinedsex, RNA_Lip_auto_combinedsex, affy_Lip_auto_combinedsex)
Overlap_Chart_Lip(DNAm_Lip_auto_F, RNA_Lip_auto_F, affy_Lip_auto_F)
Overlap_Chart_Lip(DNAm_Lip_auto_M, RNA_Lip_auto_M, affy_Lip_auto_M)
Overlap_Chart_Lip(DNAm_Lip_chrX_F, RNA_Lip_chrX_F, affy_Lip_chrX_F)
Overlap_Chart_Lip(DNAm_Lip_chrX_M, RNA_Lip_chrX_M, affy_Lip_chrX_M)

Overlap_Chart_WG <- function(DNAm_file, RNA_seq_file, affy_file){
    affy_name <- deparse(substitute(affy_file))
    no_affy_name <- gsub("affy_", "", affy_name)
    filename <- paste0("overlap_", no_affy_name, ".csv")
    print(filename)
    
    DNAm <- DNAm_file[, colnames(DNAm_file) %in% c("gene", "Closest_TSS_gene_name", "probes", "gene_id", "logFC", "adj.P.Val", "deltaB", "region_overlap", "chr")]
    sig_DNAm <- DNAm[DNAm$adj.P.Val < 0.05,]
    RNA_Seq <- RNA_seq_file[, colnames(RNA_seq_file) %in% c("Row.names","gene_name", "log2FoldChange", "padj", "Expression_Status")]
    sig_RNA_Seq <- RNA_Seq[RNA_Seq$padj < 0.05,]
    affy <- affy_file[, colnames(affy_file) %in% c("ENSEMBL", "SYMBOL", "adj.P.Val", "logFC", "deltaExprs")]

    RNAmerge <- merge(sig_RNA_Seq, affy, by.x = "gene_name", by.y = "SYMBOL", all = TRUE)
    overlap_output <- merge(RNAmerge, DNAm, by.x = "gene_name", by.y = "gene", all = TRUE)
    overlap_output <- merge(overlap_output, gtf_df[ , c("gene_id", "gene_type")], by = "gene_id")

    write.csv(overlap_output, file = filename)
}

DNAm_WG_auto_combinedsex <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_DNAm_Lipid_Candidate/WG_Diff_DNAm_Analysis_2025/wg_placmet_wholepop_auto_rerun.csv")
DNAm_WG_auto_F <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_DNAm_Lipid_Candidate/WG_Diff_DNAm_Analysis_2025/wg_placmet_F_fulldata_auto_rerun.csv")
DNAm_WG_auto_M <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_DNAm_Lipid_Candidate/WG_Diff_DNAm_Analysis_2025/wg_placmet_M_fulldata_auto_rerun.csv")
DNAm_WG_chrX_F <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_DNAm_Lipid_Candidate/WG_Diff_DNAm_Analysis_2025/wg_placmet_F_fulldata_X_rerun.csv")
DNAm_WG_chrX_M <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_DNAm_Lipid_Candidate/WG_Diff_DNAm_Analysis_2025/wg_placmet_M_fulldata_X_rerun.csv")
DNAm_WG_chrY_M <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_DNAm_Lipid_Candidate/WG_Diff_DNAm_Analysis_2025/wg_placmet_M_fulldata_Y_rerun.csv")

RNA_WG_auto_combinedsex <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_WG_rerun/RNA_WG_DESeq_results_autosomes_combined_sex.csv")
RNA_WG_auto_F <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_WG_rerun/RNA_WG_DESeq_results_autosomes_F.csv")
RNA_WG_auto_M <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_WG_rerun/RNA_WG_DESeq_results_autosomes_M.csv")
RNA_WG_chrX_F <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_WG_rerun/RNA_WG_DESeq_results_chrX_F.csv")
RNA_WG_chrX_M <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_WG_rerun/RNA_WG_DESeq_results_chrX_M.csv")
RNA_WG_chrY_M <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_WG_rerun/RNA_WG_DESeq_results_chrY_M.csv")

affy_WG_auto_combinedsex <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_WG_rerun/affymetrix_validation_WG_rerun/affy_results_auto_combined_sex.csv") 
affy_WG_auto_F <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_WG_rerun/affymetrix_validation_WG_rerun/affy_results_auto_F.csv")
affy_WG_auto_M <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_WG_rerun/affymetrix_validation_WG_rerun/affy_results_auto_M.csv")
affy_WG_chrX_F <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_WG_rerun/affymetrix_validation_WG_rerun/affy_results_chrX_F.csv")
affy_WG_chrX_M <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_WG_rerun/affymetrix_validation_WG_rerun/affy_results_chrX_M.csv")

Overlap_Chart_WG(DNAm_WG_auto_combinedsex, RNA_WG_auto_combinedsex, affy_WG_auto_combinedsex)
Overlap_Chart_WG(DNAm_WG_auto_F, RNA_WG_auto_F, affy_WG_auto_F)
Overlap_Chart_WG(DNAm_WG_auto_M, RNA_WG_auto_M, affy_WG_auto_M)
Overlap_Chart_WG(DNAm_WG_chrX_F, RNA_WG_chrX_F, affy_WG_chrX_F)
Overlap_Chart_WG(DNAm_WG_chrX_M, RNA_WG_chrX_M, affy_WG_chrX_M)



# DNAm

# c("gene", "Closest_TSS_gene_name", "probe", "gene_id", "logFC", "adj.P.Val", "deltaB", "region_overlap")

# RNA-Seq
# #Row.names = gene_id
# c("Row.names","gene_symbol", "log2Fold", "padj", "Expression_Status")

# Affy
# #ENSEMBL = Row.names = gene_id
# c("ENSEMBL", "SYMBOL", "adj.P.Val", "logFC", "deltaExprs")

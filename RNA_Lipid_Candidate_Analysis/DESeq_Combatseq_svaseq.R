#Step 1: Match Lipid-related gene symbols to ensembl IDs using gencode file
BiocManager::install("rtracklayer")
BiocManager::install("data.table")
BiocManager::install("DESeq2")
library(rtracklayer)
library(data.table)
library(dplyr)
library(tidyverse)
library(sva)
library(DESeq2)

setwd("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_WG_rerun/test")

gtf_file <- "../../2025_RNA_Lipid_Candidate/genome_mapping/gencode.v48.primary_assembly.annotation.gtf"
gtf_data <- import(gtf_file)
gtf_df <- as.data.frame(gtf_data) #4119244 27
gtf_df_dups_rem <- gtf_df[!duplicated(gtf_df$gene_id),] #78894 27

setDT(gtf_df_dups_rem)

#remove lipid-related genes that didn't match ensembl file, assign lipid-related genes to chromosomes
ensembl_list <- gtf_df_dups_rem[!is.na(gtf_df_dups_rem$gene_id),]
table(ensembl_list$seqnames) ##2955 genes on chrX, 672 on chrY
chrX <- ensembl_list [ensembl_list$seqnames == "chrX", ]
chrY <- ensembl_list [ensembl_list$seqnames == "chrY", ]

autosomal_chromosomes <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22")
autosomes <- ensembl_list [ensembl_list$seqnames %in% autosomal_chromosomes, ] #75022

#Step 2: load read_counts file, re-format so that column names are sample names instead of file location
read_counts <- read.table("../../2025_RNA_Lipid_Candidate/read_counts/read_counts.txt", header = TRUE)
colnames(read_counts)[7:ncol(read_counts)] <- sub("...genome_mapping.markeddup_BAMs.", "", colnames(read_counts)[7:ncol(read_counts)])
colnames(read_counts)[7:ncol(read_counts)] <- sub("_markdup.bam", "", colnames(read_counts)[7:ncol(read_counts)])
rownames(read_counts) <- read_counts$Geneid

#Load metadata sheet to only analyze non-excluded samples 
RNA_Lipid_Candidate_Metadata <- read.csv("../../2025_RNA_Lipid_Candidate/metadata/RNA_Lipid_Candidate_Metadata.csv") #442, 24
RNA_Lipid_Candidate_Metadata$exclude[RNA_Lipid_Candidate_Metadata$GSE_number == "GSE148241"] <- "exclude" #24
RNA_Metadata_ex_rem <- RNA_Lipid_Candidate_Metadata[RNA_Lipid_Candidate_Metadata$exclude != "exclude",] #122, 24

#Condense read_counts to only exclude "exclude" samples from metadata sheet, and GSE234729, and separate by chromosome
excluded_GSE <- RNA_Lipid_Candidate_Metadata[RNA_Lipid_Candidate_Metadata$GSE_number %in% c("GSE279757","GSE234729", "GSE218039", "GSE148241"), ]

filt_read_counts <- read_counts[ , !colnames(read_counts) %in% excluded_GSE$Run]
filt_read_counts <- filt_read_counts [, 7:ncol(filt_read_counts)]

# #Control for Batch Effects using ComBat-Seq
# sex_numeric <- as.numeric(as.factor(RNA_Metadata_ex_rem$predicted_fetal_sex))
# covar_mat <- cbind(sex_numeric)
# filt_read_counts_ComBat <- ComBat_seq(as.matrix(filt_read_counts), batch=RNA_Metadata_ex_rem$GSE_number, group=RNA_Metadata_ex_rem$disease_group, covar_mod=covar_mat)

# filt_read_counts_autosomes <- filt_read_counts_ComBat [rownames(filt_read_counts_ComBat) %in% autosomes$gene_id,]
# filt_read_counts_chrX <- filt_read_counts_ComBat [rownames(filt_read_counts_ComBat) %in% chrX$gene_id,]
# filt_read_counts_chrY <- filt_read_counts_ComBat [rownames(filt_read_counts_ComBat) %in% chrY$gene_id,]

# #Controlling for Batch Effects using sva-seq
# keep_sva <- rowSums(filt_read_counts > 0) > (0.7 * ncol(filt_read_counts))
# filt_counts_sva <- as.matrix(filt_read_counts[keep_sva, ])
# row_vars <- apply(filt_counts_sva, 1, var)
# filt_counts_sva <- filt_counts_sva[row_vars > 0, ]

# mod <- model.matrix(~ disease_group + predicted_fetal_sex, data = RNA_Metadata_ex_rem)
# mod0 <- model.matrix(~ predicted_fetal_sex, data = RNA_Metadata_ex_rem)
# sv_obj <- svaseq(as.matrix(filt_counts_sva), mod, mod0)

# sv_data <- as.data.frame(sv_obj$sv)
# colnames(sv_data) <- paste0("SV", 1:6)
# metadata_combined <- cbind(RNA_Metadata_ex_rem, sv_data)

filt_read_counts_autosomes <- filt_read_counts[rownames(filt_read_counts) %in% autosomes$gene_id,]
filt_read_counts_chrX <- filt_read_counts[rownames(filt_read_counts) %in% chrX$gene_id,]
filt_read_counts_chrY <- filt_read_counts[rownames(filt_read_counts) %in% chrY$gene_id,]

#Rownames in metadata match with colnames of filt_read_counts, and in same order
rownames(RNA_Metadata_ex_rem) <- RNA_Metadata_ex_rem$Run
all(colnames(filt_read_counts) %in% rownames(RNA_Metadata_ex_rem)) #TRUE
all(colnames(filt_read_counts) == rownames(RNA_Metadata_ex_rem)) #TRUE

#Step 3: Construct DESeqDataSet object
#RNA_WG_DESeq_auto <- DESeqDataSetFromMatrix(countData = filt_read_counts_autosomes, colData = metadata_combined, design = ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + disease_group + predicted_fetal_sex)
RNA_WG_DESeq_auto <- DESeqDataSetFromMatrix(countData = filt_read_counts_autosomes, colData = RNA_Metadata_ex_rem, design = ~ GSE_number + predicted_fetal_sex + disease_group)

#Pre-filtering: removing rows with low gene counts, keeping rows that have at least 10 reads adding across all samples
pre_filter <- rowSums(counts(RNA_WG_DESeq_auto)) >=10
RNA_WG_DESeq_auto <- RNA_WG_DESeq_auto[pre_filter,]
#RNA_Lip_DESeq_auto <- RNA_WG_DESeq_auto[rownames(RNA_WG_DESeq_auto) %in% lipid_ensembl_list$gene_id]

#Set the factor level
RNA_WG_DESeq_auto$disease_group <- relevel(RNA_WG_DESeq_auto$disease_group, ref = "Control")

#Step 4: Run DESeq
RNA_WG_DESeq_auto <- DESeq(RNA_WG_DESeq_auto)
RNA_WG_DESeq_auto_results <- results(RNA_WG_DESeq_auto)
RNA_WG_DESeq_auto_results <- as.data.frame(RNA_WG_DESeq_auto_results)
RNA_WG_DESeq_auto_results <- merge(RNA_WG_DESeq_auto_results, ensembl_list[, c("gene_id", "gene_name", "seqnames")], by.x = "row.names", by.y = "gene_id")
RNA_WG_DESeq_auto_results <- RNA_WG_DESeq_auto_results[!is.na(RNA_WG_DESeq_auto_results$padj),]

#Labels 
RNA_WG_DESeq_auto_results$Expression_Status <- "Not_Biologically_Significant" #31935
RNA_WG_DESeq_auto_results$Expression_Status[RNA_WG_DESeq_auto_results$log2FoldChange > 0.00 & RNA_WG_DESeq_auto_results$padj <0.05] <- "Trending_Towards_Increased_RNA_Expression" #24
RNA_WG_DESeq_auto_results$Expression_Status[RNA_WG_DESeq_auto_results$log2FoldChange < 0.00 & RNA_WG_DESeq_auto_results$padj <0.05] <- "Trending_Towards_Decreased_RNA_Expression" #16
RNA_WG_DESeq_auto_results$Expression_Status[RNA_WG_DESeq_auto_results$log2FoldChange > 1.00 & RNA_WG_DESeq_auto_results$padj <0.05] <- "Increased_RNA_Expression" #62
RNA_WG_DESeq_auto_results$Expression_Status[RNA_WG_DESeq_auto_results$log2FoldChange < -1.00 & RNA_WG_DESeq_auto_results$padj <0.05] <- "Decreased_RNA_Expression" #30

write.csv(RNA_WG_DESeq_auto_results, "RNA_WG_DESeq_results_autosomes_combined_sex.csv")


#Female-Specific Differential Expression Analysis, autosomes
Female_Samples <- RNA_Metadata_ex_rem[RNA_Metadata_ex_rem$predicted_fetal_sex == "F",]
rownames(Female_Samples) <- Female_Samples$Run
filt_read_counts_F<- filt_read_counts[,colnames(filt_read_counts) %in% rownames(Female_Samples)]

# keep_sva_F <- rowSums(filt_read_counts > 0) > (0.7 * ncol(filt_read_counts_F))
# filt_counts_sva_F <- as.matrix(filt_read_counts_F[keep_sva_F, ])
# row_vars_F <- apply(filt_counts_sva_F, 1, var)
# filt_counts_sva_F <- filt_counts_sva_F[row_vars_F > 0, ]

# mod_F <- model.matrix(~ disease_group, data = Female_Samples)
# mod0_F <- model.matrix(~ 1, data = Female_Samples)
# sv_obj_F <- svaseq(as.matrix(filt_counts_sva_F), mod_F, mod0_F)

# sv_data_F <- as.data.frame(sv_obj_F$sv)
# colnames(sv_data_F) <- paste0("SV", 1:3)
# metadata_F <- cbind(Female_Samples, sv_data_F)


# filt_read_counts_ComBat_F <- ComBat_seq(as.matrix(filt_read_counts_F), batch=Female_Samples$GSE_number, group=Female_Samples$disease_group)
filt_read_counts_auto_F <- filt_read_counts_F[rownames(filt_read_counts_F) %in% autosomes$gene_id,]
filt_read_counts_chrX_F <- filt_read_counts_F[rownames(filt_read_counts_F) %in% chrX$gene_id,]

all(colnames(filt_read_counts_auto_F) %in% rownames(Female_Samples)) #TRUE
all(colnames(filt_read_counts_auto_F) == rownames(Female_Samples)) #TRUE

RNA_WG_DESeq_auto_F <- DESeqDataSetFromMatrix(countData = filt_read_counts_auto_F, colData = Female_Samples, design = ~ GSE_number + disease_group)

pre_filter <- rowSums(counts(RNA_WG_DESeq_auto_F)) >=10
RNA_WG_DESeq_auto_F <- RNA_WG_DESeq_auto_F[pre_filter,]
#RNA_Lip_DESeq_auto_F <- RNA_WG_DESeq_auto_F[rownames(RNA_WG_DESeq_auto_F) %in% lipid_ensembl_list$gene_id]

RNA_WG_DESeq_auto_F$disease_group <- relevel(RNA_WG_DESeq_auto_F$disease_group, ref = "Control")
RNA_WG_DESeq_auto_F <- DESeq(RNA_WG_DESeq_auto_F)

RNA_WG_DESeq_auto_F_results <- results(RNA_WG_DESeq_auto_F)
RNA_WG_DESeq_auto_F_results <- as.data.frame(RNA_WG_DESeq_auto_F_results)
RNA_WG_DESeq_auto_F_results <- merge(RNA_WG_DESeq_auto_F_results, ensembl_list[, c("gene_id", "gene_name", "seqnames")], by.x = "row.names", by.y = "gene_id")
RNA_WG_DESeq_auto_F_results <- RNA_WG_DESeq_auto_F_results[!is.na(RNA_WG_DESeq_auto_F_results$padj),]

#Labels 
RNA_WG_DESeq_auto_F_results$Expression_Status <- "Not_Biologically_Significant" #31728
RNA_WG_DESeq_auto_F_results$Expression_Status[RNA_WG_DESeq_auto_F_results$log2FoldChange > 0.00 & RNA_WG_DESeq_auto_F_results$padj <0.05] <- "Trending_Towards_Increased_RNA_Expression" #0
RNA_WG_DESeq_auto_F_results$Expression_Status[RNA_WG_DESeq_auto_F_results$log2FoldChange < 0.00 & RNA_WG_DESeq_auto_F_results$padj <0.05] <- "Trending_Towards_Decreased_RNA_Expression" #0
RNA_WG_DESeq_auto_F_results$Expression_Status[RNA_WG_DESeq_auto_F_results$log2FoldChange > 1.00 & RNA_WG_DESeq_auto_F_results$padj <0.05] <- "Increased_RNA_Expression" #98
RNA_WG_DESeq_auto_F_results$Expression_Status[RNA_WG_DESeq_auto_F_results$log2FoldChange < -1.00 & RNA_WG_DESeq_auto_F_results$padj <0.05] <- "Decreased_RNA_Expression" #176

write.csv(RNA_WG_DESeq_auto_F_results, "RNA_WG_DESeq_results_autosomes_F.csv")


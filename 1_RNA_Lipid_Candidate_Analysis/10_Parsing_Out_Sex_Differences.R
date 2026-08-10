#Identify which differentially methylated genes are common and unique to male and female

##Lipid-Candidate

affy_Lip_auto_F <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/affymetrix_validation_rerun/affy_results_auto_F.csv")
affy_Lip_auto_M <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/affymetrix_validation_rerun/affy_results_auto_M.csv")
bio_sig_affy_F <- affy_Lip_auto_F[affy_Lip_auto_F$adj.P.Val <= 0.05 & (affy_Lip_auto_F$logFC >= 1 | affy_Lip_auto_F$logFC <= -1), ]
sig_affy_F <- affy_Lip_auto_F[affy_Lip_auto_F$adj.P.Val <= 0.05, ]
bio_sig_affy_M <- affy_Lip_auto_M[affy_Lip_auto_M$adj.P.Val <= 0.05 & (affy_Lip_auto_M$logFC >= 1 | affy_Lip_auto_M$logFC <= -1), ]
sig_affy_M <- affy_Lip_auto_M[affy_Lip_auto_M$adj.P.Val <= 0.05, ]

#Genes sig only in females and only iin males - affy
only_bio_sig_affy_F <- sig_affy_F[!(sig_affy_F$SYMBOL %in% sig_affy_M$SYMBOL),] #158 genes
only_bio_sig_affy_M <- sig_affy_M[!(sig_affy_M$SYMBOL %in% sig_affy_F$SYMBOL),] #96 genes
#genes in both but with pvals + logFC in females, then males - affy
female_male_DEG_affy_F <- sig_affy_F[(sig_affy_F$SYMBOL %in% sig_affy_M$SYMBOL),] #259 genes
female_male_DEG_affy_M <- sig_affy_M[(sig_affy_M$SYMBOL %in% sig_affy_F$SYMBOL),]

affy_Lip_auto_F$unique <- "Neither"
affy_Lip_auto_F$unique[affy_Lip_auto_F$SYMBOL %in% only_bio_sig_affy_F$SYMBOL] <- "Only_Female"
affy_Lip_auto_F$unique[affy_Lip_auto_F$SYMBOL %in% female_male_DEG_affy_F$SYMBOL] <- "Both"

affy_Lip_auto_M$unique <- "Neither"
affy_Lip_auto_M$unique[affy_Lip_auto_M$SYMBOL %in% only_bio_sig_affy_M$SYMBOL] <- "Only_Male"
affy_Lip_auto_M$unique[affy_Lip_auto_M$SYMBOL %in% female_male_DEG_affy_M$SYMBOL] <- "Both"


RNA_Lip_auto_F <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/DESeq_rerun/RNA_Lip_DESeq_results_autosomes_F.csv")
RNA_Lip_auto_M <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/DESeq_rerun/RNA_Lip_DESeq_results_autosomes_M.csv")
bio_sig_F <- RNA_Lip_auto_F[RNA_Lip_auto_F$Expression_Status %in% (c("Decreased_RNA_Expression", "Increased_RNA_Expression")), ]
bio_sig_M <- RNA_Lip_auto_M[RNA_Lip_auto_M$Expression_Status %in% (c("Decreased_RNA_Expression", "Increased_RNA_Expression")), ]

#Genes sig only in females and only iin males - RNAseq
only_bio_sig_F <- bio_sig_F[!(bio_sig_F$gene_symbol %in% bio_sig_M$gene_symbol),] #173 genes
only_bio_sig_M <- bio_sig_M[!(bio_sig_M$gene_symbol %in% bio_sig_F$gene_symbol),] #242 genes
#genes in both but with pvals + logFC in females, then males - RNAseq
female_male_DEG_F <- bio_sig_F[(bio_sig_F$gene_symbol %in% bio_sig_M$gene_symbol),] #259 genes
female_male_DEG_M <- bio_sig_M[(bio_sig_M$gene_symbol %in% bio_sig_F$gene_symbol),]

RNA_Lip_auto_F$unique <- "Neither"
RNA_Lip_auto_F$unique[RNA_Lip_auto_F$gene_symbol %in% only_bio_sig_F$gene_symbol] <- "Only_Female"
RNA_Lip_auto_F$unique[RNA_Lip_auto_F$gene_symbol %in% female_male_DEG_F$gene_symbol] <- "Both"

RNA_Lip_auto_M$unique <- "Neither"
RNA_Lip_auto_M$unique[RNA_Lip_auto_M$gene_symbol %in% only_bio_sig_M$gene_symbol] <- "Only_Male"
RNA_Lip_auto_M$unique[RNA_Lip_auto_M$gene_symbol %in% female_male_DEG_M$gene_symbol] <- "Both"

RNAmerge_Lip_onlyF <- merge(RNA_Lip_auto_F, affy_Lip_auto_F, by.x = "gene_symbol", by.y = "SYMBOL", all = TRUE)
RNAmerge_Lip_onlyM <- merge(RNA_Lip_auto_M, affy_Lip_auto_M, by.x = "gene_symbol", by.y = "SYMBOL", all = TRUE)


write.csv(RNA_Lip_auto_F, "RNA_Lip_DESeq_results_autosomes_F.csv")
write.csv(RNA_Lip_auto_M, "RNA_Lip_DESeq_results_autosomes_M.csv")

write.csv(RNAmerge_Lip_onlyF, "RNAmerge_Lip_autosomes_onlyF.csv")
write.csv(RNAmerge_Lip_onlyM, "RNAmerge_Lip_autosomes_onlyM.csv")

##Whole Genome
affy_WG_auto_F <-read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_WG_rerun/affymetrix_validation_WG_rerun/affy_results_auto_F.csv")
affy_WG_auto_M <-read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_WG_rerun/affymetrix_validation_WG_rerun/affy_results_auto_M.csv")

RNA_Lip_auto_F_WG <-read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_WG_rerun/RNA_WG_DESeq_results_autosomes_F.csv")
RNA_Lip_auto_M_WG <-read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_WG_rerun/RNA_WG_DESeq_results_autosomes_M.csv")
bio_sig_F <- RNA_Lip_auto_F_WG[RNA_Lip_auto_F_WG$Expression_Status %in% (c("Decreased_RNA_Expression", "Increased_RNA_Expression")), ]
bio_sig_M <- RNA_Lip_auto_M_WG[RNA_Lip_auto_M_WG$Expression_Status %in% (c("Decreased_RNA_Expression", "Increased_RNA_Expression")), ]
bio_sig_affy_F <- affy_WG_auto_F[affy_WG_auto_F$adj.P.Val <= 0.05 & (affy_WG_auto_F$logFC >= 1 | affy_WG_auto_F$logFC <= -1), ]
sig_affy_F <- affy_WG_auto_F[affy_WG_auto_F$adj.P.Val <= 0.05, ]
bio_sig_affy_M <- affy_WG_auto_M[affy_WG_auto_M$adj.P.Val <= 0.05 & (affy_WG_auto_M$logFC >= 1 | affy_WG_auto_M$logFC <= -1), ]
sig_affy_M <- affy_WG_auto_M[affy_WG_auto_M$adj.P.Val <= 0.05, ]

#Genes sig only in females and only iin males - affy
only_bio_sig_affy_F <- sig_affy_F[!(sig_affy_F$SYMBOL %in% sig_affy_M$SYMBOL),] #158 genes
only_bio_sig_affy_M <- sig_affy_M[!(sig_affy_M$SYMBOL %in% sig_affy_F$SYMBOL),] #96 genes
#genes in both but with pvals + logFC in females, then males - affy
female_male_DEG_affy_F <- sig_affy_F[(sig_affy_F$SYMBOL %in% sig_affy_M$SYMBOL),] #259 genes
female_male_DEG_affy_M <- sig_affy_M[(sig_affy_M$SYMBOL %in% sig_affy_F$SYMBOL),]

affy_WG_auto_F$unique <- "Neither"
affy_WG_auto_F$unique[affy_WG_auto_F$SYMBOL %in% only_bio_sig_affy_F$SYMBOL] <- "Only_Female"
affy_WG_auto_F$unique[affy_WG_auto_F$SYMBOL %in% female_male_DEG_affy_F$SYMBOL] <- "Both"

affy_WG_auto_M$unique <- "Neither"
affy_WG_auto_M$unique[affy_WG_auto_M$SYMBOL %in% only_bio_sig_affy_M$SYMBOL] <- "Only_Male"
affy_WG_auto_M$unique[affy_WG_auto_M$SYMBOL %in% female_male_DEG_affy_M$SYMBOL] <- "Both"

#Genes sig only in females and only iin males
only_bio_sig_F <- bio_sig_F[!(bio_sig_F$gene_name %in% bio_sig_M$gene_name),] 
only_bio_sig_M <- bio_sig_M[!(bio_sig_M$gene_name %in% bio_sig_F$gene_name),] 
#genes in both but with pvals + logFC in females, then males
female_male_DEG_F <- bio_sig_F[(bio_sig_F$gene_name %in% bio_sig_M$gene_name),] 
female_male_DEG_M <- bio_sig_M[(bio_sig_M$gene_name %in% bio_sig_F$gene_name),]

RNA_Lip_auto_F_WG$unique <- "Neither"
RNA_Lip_auto_F_WG$unique[RNA_Lip_auto_F_WG$gene_name %in% only_bio_sig_F$gene_name] <- "Only_Female"
RNA_Lip_auto_F_WG$unique[RNA_Lip_auto_F_WG$gene_name %in% female_male_DEG_F$gene_name] <- "Both"

RNA_Lip_auto_M_WG$unique <- "Neither"
RNA_Lip_auto_M_WG$unique[RNA_Lip_auto_M_WG$gene_name %in% only_bio_sig_M$gene_name] <- "Only_Male"
RNA_Lip_auto_M_WG$unique[RNA_Lip_auto_M_WG$gene_name %in% female_male_DEG_M$gene_name] <- "Both"

write.csv(RNA_Lip_auto_F_WG, "RNA_Lip_DESeq_results_autosomes_F.csv")
write.csv(RNA_Lip_auto_M_WG, "RNA_Lip_DESeq_results_autosomes_M.csv")

RNAmerge_WG_onlyF <- merge(RNA_Lip_auto_F_WG, affy_WG_auto_F, by.x = "gene_name", by.y = "SYMBOL", all = TRUE)
RNAmerge_WG_onlyM <- merge(RNA_Lip_auto_M_WG, affy_WG_auto_M, by.x = "gene_name", by.y = "SYMBOL", all = TRUE)

write.csv(RNAmerge_WG_onlyF, "RNAmerge_WG_autosomes_onlyF.csv")
write.csv(RNAmerge_WG_onlyM, "RNAmerge_WG_autosomes_onlyM.csv")

#Rerun DESeq2 code (combined_sex) to include interaction term between fetal-sex and disease group

#Step 1: Match Lipid-related gene symbols to ensembl IDs using gencode file
BiocManager::install("rtracklayer")
BiocManager::install("data.table")
BiocManager::install("DESeq2")
library(rtracklayer)
library(data.table)
library(dplyr)
library(tidyverse)
library(DESeq2)

setwd("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/DESeq_rerun")

gtf19_file <- "./gencode.v37lift37.annotation.gtf"
gtf19_data <- import(gtf19_file)
gtf19_df <- as.data.frame(gtf19_data)

gtf_file <- "../genome_mapping/gencode.v48.primary_assembly.annotation.gtf"
gtf_data <- import(gtf_file)
gtf_df <- as.data.frame(gtf_data)
#gtf_df$gene_id <- sub("\\..*", "", gtf_df$gene_id)
gtf_df_dups_rem <- gtf_df[!duplicated(gtf_df$gene_id),]

lipid_gene_list <- read.csv("lipid_genes_unique_P.csv")

setDT(lipid_gene_list)
setDT(gtf_df_dups_rem)
lipid_ensembl_gencode <- merge(x = lipid_gene_list, y = gtf_df_dups_rem[, c("gene_id", "gene_name", "seqnames")], by.x = "gene_symbol", by.y = "gene_name", all.x = TRUE)
#lipid_ensembl <- named_read_counts[,c("Geneid", "gene_name")]
write.csv(lipid_ensembl_gencode, "lipid_ensembl_gencode.csv")

#remove lipid-related genes that didn't match ensembl file, assign lipid-related genes to chromosomes
lipid_ensembl_list <- lipid_ensembl_gencode[!is.na(lipid_ensembl_gencode$gene_id),]
table(lipid_ensembl_list$seqnames) ##170 lipid-related genes on chrX, 7 on chrY
chrX_lipids <- lipid_ensembl_list [lipid_ensembl_list$seqnames == "chrX", ]
chrY_lipids <- lipid_ensembl_list [lipid_ensembl_list$seqnames == "chrY", ]

autosomal_chromosomes <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22")
autosomal_lipids <- lipid_ensembl_list [lipid_ensembl_list$seqnames %in% autosomal_chromosomes, ]

#Step 2: load read_counts file, re-format so that column names are sample names instead of file location
read_counts <- read.table("../read_counts/read_counts.txt", header = TRUE)
colnames(read_counts)[7:ncol(read_counts)] <- sub("...genome_mapping.markeddup_BAMs.", "", colnames(read_counts)[7:ncol(read_counts)])
colnames(read_counts)[7:ncol(read_counts)] <- sub("_markdup.bam", "", colnames(read_counts)[7:ncol(read_counts)])
rownames(read_counts) <- read_counts$Geneid

#Load metadata sheet to only analyze non-excluded samples 
RNA_Lipid_Candidate_Metadata <- read.csv("../metadata/RNA_Lipid_Candidate_Metadata.csv")
RNA_Metadata_ex_rem <- RNA_Lipid_Candidate_Metadata[RNA_Lipid_Candidate_Metadata$exclude != "exclude",]

#Condense read_counts to only exclude "exclude" samples from metadata sheet, and GSE234729, and separate by chromosome
excluded_GSE <- RNA_Lipid_Candidate_Metadata[RNA_Lipid_Candidate_Metadata$GSE_number %in% c("GSE279757","GSE234729","GSE218039"), ]

filt_read_counts <- read_counts[ , !colnames(read_counts) %in% excluded_GSE$Run]
filt_read_counts <- filt_read_counts [, 7:ncol(filt_read_counts)]

filt_read_counts_autosomes <- filt_read_counts [rownames(filt_read_counts) %in% autosomal_lipids$gene_id,]
filt_read_counts_chrX <- filt_read_counts [rownames(filt_read_counts) %in% chrX_lipids$gene_id,]
filt_read_counts_chrY <- filt_read_counts [rownames(filt_read_counts) %in% chrY_lipids$gene_id,]

#Rownames in metadata match with colnames of filt_read_counts, and in same order
rownames(RNA_Metadata_ex_rem) <- RNA_Metadata_ex_rem$Run
all(colnames(filt_read_counts) %in% rownames(RNA_Metadata_ex_rem)) #TRUE
all(colnames(filt_read_counts) == rownames(RNA_Metadata_ex_rem)) #TRUE

#Step 3: Construct DESeqDataSet object, removed "instrument" from design because confound perfectly with GSE_number
RNA_WG_DESeq_auto <- DESeqDataSetFromMatrix(countData = filt_read_counts_autosomes, colData = RNA_Metadata_ex_rem, design = ~ GSE_number + predicted_fetal_sex*disease_group)

#Pre-filtering: removing rows with low gene counts, keeping rows that have at least 10 reads adding across all samples
pre_filter <- rowSums(counts(RNA_WG_DESeq_auto)) >=10
RNA_WG_DESeq_auto <- RNA_WG_DESeq_auto[pre_filter,]
RNA_Lip_DESeq_auto <- RNA_WG_DESeq_auto[rownames(RNA_WG_DESeq_auto) %in% lipid_ensembl_list$gene_id]

#Set the factor level
RNA_Lip_DESeq_auto$disease_group <- relevel(RNA_Lip_DESeq_auto$disease_group, ref = "Control")

#Step 4: Run DESeq
RNA_Lip_DESeq_auto <- DESeq(RNA_Lip_DESeq_auto)
RNA_Lip_DESeq_auto_results <- results(RNA_Lip_DESeq_auto)
RNA_Lip_DESeq_auto_results <- as.data.frame(RNA_Lip_DESeq_auto_results)
RNA_Lip_DESeq_auto_results <- merge(RNA_Lip_DESeq_auto_results, lipid_ensembl_list[, c("gene_id", "gene_symbol", "seqnames")], by.x = "row.names", by.y = "gene_id")
RNA_Lip_DESeq_auto_results <- RNA_Lip_DESeq_auto_results[!is.na(RNA_Lip_DESeq_auto_results$padj),]

resultsNames(RNA_Lip_DESeq_auto)
res_female_disease <- results(RNA_Lip_DESeq_auto, name="disease_group_PE_vs_Control")
res_male_disease <- results(RNA_Lip_DESeq_auto, contrast=list(c("disease_group_PE_vs_Control", "predicted_fetal_sexM.disease_groupPE")))
res_interaction <- results(RNA_Lip_DESeq_auto, name="predicted_fetal_sexM.disease_groupPE")

res_female <- as.data.frame(res_female_disease)
sig_female <- res_female[!is.na(res_female$padj) & res_female$padj < 0.05, ]
sig_female$comparison <- "F_PEvsF_Cont" # positive logFC = lipid gene expression goes up in female PE, negative logFC = lipid gene expression goes down in female PE
sig_female$gene_id <- rownames(sig_female)

all_female <- res_female[!is.na(res_female$padj),]
all_female$comparison <- "F_PEvsF_Cont" # positive logFC = lipid gene expression goes up in female PE, negative logFC = lipid gene expression goes down in female PE
all_female$gene_id <- rownames(all_female)

res_male <- as.data.frame(res_male_disease)
sig_male <- res_male[!is.na(res_male$padj) & res_male$padj < 0.05, ]
sig_male$comparison <- "M_PEvsM_Cont" # positive logFC = lipid gene expression goes up in female PE, negative logFC = lipid gene expression goes down in female PE
sig_male$gene_id <- rownames(sig_male)

all_male <- res_male[!is.na(res_male$padj), ]
all_male$comparison <- "M_PEvsM_Cont" # positive logFC = lipid gene expression goes up in female PE, negative logFC = lipid gene expression goes down in female PE
all_male$gene_id <- rownames(all_male)

sig_interaction <- as.data.frame(res_interaction)
sig_interaction <- sig_interaction[!is.na(sig_interaction$padj) & sig_interaction$padj < 0.05, ]
sig_interaction$comparison <- "Interaction_PEvsF_Cont" # positive logFC = more of increase or less of decrease in gene expression in Male PE than Female PE, negative logFC = more of gene decrease, less of gene expression increase in Male PE than Female PE
sig_interaction$gene_id <- rownames(sig_interaction)

Lip_IT_result <- bind_rows(sig_female,sig_male,sig_interaction)
Lip_IT_result <- merge(Lip_IT_result, gtf_df_dups_rem[, c("gene_id", "gene_name")], by = "gene_id", all.x = TRUE)

SS_IT_result <- bind_rows(all_female,all_male)
SS_IT_result <- merge(SS_IT_result, gtf_df_dups_rem[, c("gene_id", "gene_name")], by = "gene_id", all.x = TRUE)

#Labels 
Lip_IT_result$Expression_Status <- "Not_Biologically_Significant" #4753
Lip_IT_result$Expression_Status[Lip_IT_result$log2FoldChange > 0.00 & Lip_IT_result$padj <0.05] <- "Trending_Towards_Increased_RNA_Expression" #0
Lip_IT_result$Expression_Status[Lip_IT_result$log2FoldChange < 0.00 & Lip_IT_result$padj <0.05] <- "Trending_Towards_Decreased_RNA_Expression" #0
Lip_IT_result$Expression_Status[Lip_IT_result$log2FoldChange > 1.00 & Lip_IT_result$padj <0.05] <- "Increased_RNA_Expression" #7
Lip_IT_result$Expression_Status[Lip_IT_result$log2FoldChange < -1.00 & Lip_IT_result$padj <0.05] <- "Decreased_RNA_Expression" #7

SS_IT_result$Expression_Status <- "Not_Biologically_Significant" #4753
SS_IT_result$Expression_Status[SS_IT_result$log2FoldChange > 0.00 & SS_IT_result$padj <0.05] <- "Trending_Towards_Increased_RNA_Expression" #0
SS_IT_result$Expression_Status[SS_IT_result$log2FoldChange < 0.00 & SS_IT_result$padj <0.05] <- "Trending_Towards_Decreased_RNA_Expression" #0
SS_IT_result$Expression_Status[SS_IT_result$log2FoldChange > 1.00 & SS_IT_result$padj <0.05] <- "Increased_RNA_Expression" #7
SS_IT_result$Expression_Status[SS_IT_result$log2FoldChange < -1.00 & SS_IT_result$padj <0.05] <- "Decreased_RNA_Expression" #7

write.csv(Lip_IT_result, "RNA_Lip_DESeq_results_autosomes_IT.csv")
write.csv(SS_IT_result, "RNA_Lip_DESeq_results_autosomes_IT_all.csv")


##Whole-Genome
setwd("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_WG_rerun")

#remove lipid-related genes that didn't match ensembl file, assign lipid-related genes to chromosomes
ensembl_list <- gtf_df_dups_rem[!is.na(gtf_df_dups_rem$gene_id),]
table(ensembl_list$seqnames) ##2955 genes on chrX, 672 on chrY
chrX <- ensembl_list [ensembl_list$seqnames == "chrX", ]
chrY <- ensembl_list [ensembl_list$seqnames == "chrY", ]

autosomal_chromosomes <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22")
autosomes <- ensembl_list [ensembl_list$seqnames %in% autosomal_chromosomes, ] #75022

#Step 2: load read_counts file, re-format so that column names are sample names instead of file location
read_counts <- read.table("../2025_RNA_Lipid_Candidate/read_counts/read_counts.txt", header = TRUE)
colnames(read_counts)[7:ncol(read_counts)] <- sub("...genome_mapping.markeddup_BAMs.", "", colnames(read_counts)[7:ncol(read_counts)])
colnames(read_counts)[7:ncol(read_counts)] <- sub("_markdup.bam", "", colnames(read_counts)[7:ncol(read_counts)])
rownames(read_counts) <- read_counts$Geneid

#Load metadata sheet to only analyze non-excluded samples 
RNA_Lipid_Candidate_Metadata <- read.csv("../2025_RNA_Lipid_Candidate/metadata/RNA_Lipid_Candidate_Metadata.csv") #442, 24
RNA_Metadata_ex_rem <- RNA_Lipid_Candidate_Metadata[RNA_Lipid_Candidate_Metadata$exclude != "exclude",] #122, 24

#Condense read_counts to only exclude "exclude" samples from metadata sheet, and GSE234729, and separate by chromosome
excluded_GSE <- RNA_Lipid_Candidate_Metadata[RNA_Lipid_Candidate_Metadata$GSE_number %in% c("GSE279757","GSE234729", "GSE218039"), ]

filt_read_counts <- read_counts[ , !colnames(read_counts) %in% excluded_GSE$Run]
filt_read_counts <- filt_read_counts [, 7:ncol(filt_read_counts)]

filt_read_counts_autosomes <- filt_read_counts[rownames(filt_read_counts) %in% autosomes$gene_id,]
filt_read_counts_chrX <- filt_read_counts[rownames(filt_read_counts) %in% chrX$gene_id,]
filt_read_counts_chrY <- filt_read_counts[rownames(filt_read_counts) %in% chrY$gene_id,]

#Rownames in metadata match with colnames of filt_read_counts, and in same order
rownames(RNA_Metadata_ex_rem) <- RNA_Metadata_ex_rem$Run
all(colnames(filt_read_counts) %in% rownames(RNA_Metadata_ex_rem)) #TRUE
all(colnames(filt_read_counts) == rownames(RNA_Metadata_ex_rem)) #TRUE

#Step 3: Construct DESeqDataSet object, removed "instrument" from design because confound perfectly with GSE_number
RNA_WG_DESeq_auto <- DESeqDataSetFromMatrix(countData = filt_read_counts_autosomes, colData = RNA_Metadata_ex_rem, design = ~ GSE_number + predicted_fetal_sex*disease_group)

#Pre-filtering: removing rows with low gene counts, keeping rows that have at least 10 reads adding across all samples
pre_filter <- rowSums(counts(RNA_WG_DESeq_auto)) >=10
RNA_WG_DESeq_auto <- RNA_WG_DESeq_auto[pre_filter,]

#Set the factor level
RNA_WG_DESeq_auto$disease_group <- relevel(RNA_WG_DESeq_auto$disease_group, ref = "Control")

#Step 4: Run DESeq
RNA_WG_DESeq_auto <- DESeq(RNA_WG_DESeq_auto)

resultsNames(RNA_WG_DESeq_auto)
res_female_disease <- results(RNA_WG_DESeq_auto, name="disease_group_PE_vs_Control")
res_male_disease <- results(RNA_WG_DESeq_auto, contrast=list(c("disease_group_PE_vs_Control", "predicted_fetal_sexM.disease_groupPE")))
res_interaction <- results(RNA_WG_DESeq_auto, name="predicted_fetal_sexM.disease_groupPE")


res_female <- as.data.frame(res_female_disease)
sig_female <- res_female[!is.na(res_female$padj) & res_female$padj < 0.05, ]
sig_female$comparison <- "F_PEvsF_Cont" # positive logFC = lipid gene expression goes up in female PE, negative logFC = lipid gene expression goes down in female PE
sig_female$gene_id <- rownames(sig_female)

all_female <- res_female[!is.na(res_female$padj),]
all_female$comparison <- "F_PEvsF_Cont" # positive logFC = lipid gene expression goes up in female PE, negative logFC = lipid gene expression goes down in female PE
all_female$gene_id <- rownames(all_female)

res_male <- as.data.frame(res_male_disease)
sig_male <- res_male[!is.na(res_male$padj) & res_male$padj < 0.05, ]
sig_male$comparison <- "M_PEvsM_Cont" # positive logFC = lipid gene expression goes up in female PE, negative logFC = lipid gene expression goes down in female PE
sig_male$gene_id <- rownames(sig_male)

all_male <- res_male[!is.na(res_male$padj), ]
all_male$comparison <- "M_PEvsM_Cont" # positive logFC = lipid gene expression goes up in female PE, negative logFC = lipid gene expression goes down in female PE
all_male$gene_id <- rownames(all_male)

sig_interaction <- as.data.frame(res_interaction)
sig_interaction <- sig_interaction[!is.na(sig_interaction$padj) & sig_interaction$padj < 0.05, ]
sig_interaction$comparison <- "Interaction_PEvsF_Cont" # positive logFC = more of increase or less of decrease in gene expression in Male PE than Female PE, negative logFC = more of gene decrease, less of gene expression increase in Male PE than Female PE
sig_interaction$gene_id <- rownames(sig_interaction)

WG_IT_result <- bind_rows(sig_female,sig_male,sig_interaction)
WG_IT_result <- merge(WG_IT_result, gtf_df_dups_rem[, c("gene_id", "gene_name")], by = "gene_id", all.x = TRUE)

SS_WG_IT_result <- bind_rows(all_female,all_male)
SS_WG_IT_result <- merge(SS_WG_IT_result, gtf_df_dups_rem[, c("gene_id", "gene_name")], by = "gene_id", all.x = TRUE)

#Labels 
WG_IT_result$Expression_Status <- "Not_Biologically_Significant" #4753
WG_IT_result$Expression_Status[WG_IT_result$log2FoldChange > 0.00 & WG_IT_result$padj <0.05] <- "Trending_Towards_Increased_RNA_Expression" #0
WG_IT_result$Expression_Status[WG_IT_result$log2FoldChange < 0.00 & WG_IT_result$padj <0.05] <- "Trending_Towards_Decreased_RNA_Expression" #0
WG_IT_result$Expression_Status[WG_IT_result$log2FoldChange > 1.00 & WG_IT_result$padj <0.05] <- "Increased_RNA_Expression" #7
WG_IT_result$Expression_Status[WG_IT_result$log2FoldChange < -1.00 & WG_IT_result$padj <0.05] <- "Decreased_RNA_Expression" #7

SS_WG_IT_result$Expression_Status <- "Not_Biologically_Significant" #4753
SS_WG_IT_result$Expression_Status[SS_WG_IT_result$log2FoldChange > 0.00 & SS_WG_IT_result$padj <0.05] <- "Trending_Towards_Increased_RNA_Expression" #0
SS_WG_IT_result$Expression_Status[SS_WG_IT_result$log2FoldChange < 0.00 & SS_WG_IT_result$padj <0.05] <- "Trending_Towards_Decreased_RNA_Expression" #0
SS_WG_IT_result$Expression_Status[SS_WG_IT_result$log2FoldChange > 1.00 & SS_WG_IT_result$padj <0.05] <- "Increased_RNA_Expression" #7
SS_WG_IT_result$Expression_Status[SS_WG_IT_result$log2FoldChange < -1.00 & SS_WG_IT_result$padj <0.05] <- "Decreased_RNA_Expression" #7

write.csv(WG_IT_result, "RNA_WG_DESeq_results_autosomes_IT.csv")
write.csv(SS_WG_IT_result, "RNA_WG_DESeq_results_autosomes_IT_all.csv")


library(ggplot2)
library(dplyr)
library(gridExtra)
library(ggrepel) 

WG_IT_result <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_WG_rerun/RNA_WG_DESeq_results_autosomes_IT.csv")
WG_F_DEGs <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_WG_rerun/RNA_WG_DESeq_results_autosomes_F.csv")
colnames(WG_F_DEGs)[colnames(WG_F_DEGs) == "log2FoldChange"] <- "F_log2FoldChange"
colnames(WG_F_DEGs)[colnames(WG_F_DEGs) == "padj"] <- "F_padj"
WG_M_DEGs <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_WG_rerun/RNA_WG_DESeq_results_autosomes_M.csv")
colnames(WG_M_DEGs)[colnames(WG_M_DEGs) == "log2FoldChange"] <- "M_log2FoldChange"
colnames(WG_M_DEGs)[colnames(WG_M_DEGs) == "padj"] <- "M_padj"
WG_IT_result <- merge(WG_IT_result, WG_F_DEGs[,c("Row.names", "F_log2FoldChange", "F_padj")], by.x = "gene_id", by.y = "Row.names")
WG_IT_result <- merge(WG_IT_result, WG_M_DEGs[,c("Row.names", "M_log2FoldChange", "M_padj")], by.x = "gene_id", by.y = "Row.names")
WG_IT_result_sig <- WG_IT_result[WG_IT_result$comparison == "Interaction_PEvsF_Cont",]
WG_IT_result_sig <- WG_IT_result_sig[!WG_IT_result_sig$gene_name == "MATN3",]

F_sig <- WG_IT_result_sig$F_padj < 0.05
M_sig <- WG_IT_result_sig$M_padj < 0.05
F_up  <- WG_IT_result_sig$F_log2FoldChange > 0
M_up  <- WG_IT_result_sig$M_log2FoldChange > 0

WG_IT_result_sig$status <- "No_Significant_Sex-Differences"
WG_IT_result_sig$status[F_sig & M_sig & F_up & M_up & (abs(WG_IT_result_sig$F_log2FoldChange) > abs(WG_IT_result_sig$M_log2FoldChange))] <- "Upregulated in Both (Female > Male)"
WG_IT_result_sig$status[F_sig & M_sig & F_up & M_up & (abs(WG_IT_result_sig$F_log2FoldChange) < abs(WG_IT_result_sig$M_log2FoldChange))] <- "Upregulated in Both (Male > Female)"
WG_IT_result_sig$status[F_sig & M_sig & !F_up & !M_up & (abs(WG_IT_result_sig$F_log2FoldChange) > abs(WG_IT_result_sig$M_log2FoldChange))] <- "Downregulated in Both (Female < Male)"
WG_IT_result_sig$status[F_sig & M_sig & !F_up & !M_up & (abs(WG_IT_result_sig$F_log2FoldChange) < abs(WG_IT_result_sig$M_log2FoldChange))] <- "Downregulated in Both (Male < Female)"

WG_IT_result_sig$status[F_sig & M_sig & F_up & !M_up] <- "Opposite Expression (Female Upregulated, Male Downregulated)"
WG_IT_result_sig$status[F_sig & M_sig & !F_up & M_up] <- "Opposite Expression (Male Upregulated, Female Downregulated)"

WG_IT_result_sig$status[F_sig & !M_sig & F_up] <- "Female-Specific Upregulation (non-sig adjP in M)"
WG_IT_result_sig$status[F_sig & !M_sig & !F_up] <- "Female-Specific Downregulation (non-sig adjP in M)"

WG_IT_result_sig$status[!F_sig & M_sig & M_up] <- "Male-Specific Upregulation (non-sig adjP in F)"
WG_IT_result_sig$status[!F_sig & M_sig & !M_up] <- "Male-Specific Downregulation (non-sig adjP in F)"

# cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") +

Sex_Diff_plot <- ggplot(WG_IT_result_sig, aes(x= F_log2FoldChange, y=M_log2FoldChange, color=status)) + 
  geom_point(alpha = 0.6, size = 3) +
  scale_color_manual(values = c(
    "No_Significant_Sex-Differences" = "#000000", # Black

    "Female-Specific Upregulation (non-sig adjP in M)"  = "#E69F00", # Orange (Original Okabe-Ito)
    "Female-Specific Downregulation (non-sig adjP in M)"  = "#D55E00", # Reddish Purple / Pink (Original Okabe-Ito)
    # "Upregulated in Both (Female > Male)"   = "#D55E00", # Vermilion / Deep Red (Original Okabe-Ito)
    # "Downregulated in Both (Female < Male)"   = "#F0E442", # Yellow (Original Okabe-Ito)

    "Male-Specific Upregulation (non-sig adjP in F)"  = "#009E73", # Bluish Green (Original Okabe-Ito)
    "Male-Specific Downregulation (non-sig adjP in F)"  = "#56B4E9", # Sky Blue (Original Okabe-Ito)
    "Upregulated in Both (Male > Female)"     = "#0f62fe", # Dark Blue (Original Okabe-Ito)
    #"Downregulated in Both (Male < Female)"     = "#0f62fe", # Teal / Mint (Added Paul Tol color)

    "Opposite Expression (Female Upregulated, Male Downregulated)" = "#d02670" # Mauve / Dark Pink (Added Safe Accent)
    #"Opposite Expression (Male Upregulated, Female Downregulated)" = "#0e6027"  # Deep Forest Green (Added Safe Accent)

  
  )) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
    geom_text_repel(aes(label=gene_name), na.rm = TRUE, max.overlaps = Inf, size = 4, segment.colour = 'grey50') +
guides(color = guide_legend(nrow = 2, byrow = TRUE, title.position = "top", title.hjust = 0.5))+
  theme_minimal() +
  theme(plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 14),
        legend.position = "bottom"
        ) +
  labs(
    x = "Log2FC_Female (DEGs)",
    y = "Log2FC_Male (DEGs)",
    title = "Comparison of DEGs Between Females and Males",
    color = "Status"
  ) 

png("./M_F_DEG_Comparison.png", height = 8, width = 14, units = "in", res = 500)
print(Sex_Diff_plot)
dev.off()



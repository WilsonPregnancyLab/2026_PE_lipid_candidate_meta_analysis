
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

lipid_gene_list <- read.csv("./lipid_genes_unique_P.csv")

setDT(lipid_gene_list)
setDT(gtf_df_dups_rem)
lipid_ensembl_gencode <- merge(x = lipid_gene_list, y = gtf_df_dups_rem[, c("gene_id", "gene_name", "seqnames")], by.x = "gene_symbol", by.y = "gene_name", all.x = TRUE)
#lipid_ensembl <- named_read_counts[,c("Geneid", "gene_name")]
write.csv(lipid_ensembl_gencode, "lipid_ensembl_gencode.csv")

#remove lipid-related genes that didn't match ensembl file, assign lipid-related genes to chromosomes
lipid_ensembl_list <- lipid_ensembl_gencode[!is.na(lipid_ensembl_gencode$gene_id),] #5340 
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
RNA_WG_DESeq_auto <- DESeqDataSetFromMatrix(countData = filt_read_counts_autosomes, colData = RNA_Metadata_ex_rem, design = ~ GSE_number + predicted_fetal_sex + disease_group)

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

#Labels 
RNA_Lip_DESeq_auto_results$Expression_Status <- "Not_Biologically_Significant" #3882
RNA_Lip_DESeq_auto_results$Expression_Status[RNA_Lip_DESeq_auto_results$log2FoldChange > 0.00 & RNA_Lip_DESeq_auto_results$padj <0.05] <- "Trending_Towards_Increased_RNA_Expression" #454
RNA_Lip_DESeq_auto_results$Expression_Status[RNA_Lip_DESeq_auto_results$log2FoldChange < 0.00 & RNA_Lip_DESeq_auto_results$padj <0.05] <- "Trending_Towards_Decreased_RNA_Expression" #341
RNA_Lip_DESeq_auto_results$Expression_Status[RNA_Lip_DESeq_auto_results$log2FoldChange > 1.00 & RNA_Lip_DESeq_auto_results$padj <0.05] <- "Increased_RNA_Expression" #62
RNA_Lip_DESeq_auto_results$Expression_Status[RNA_Lip_DESeq_auto_results$log2FoldChange < -1.00 & RNA_Lip_DESeq_auto_results$padj <0.05] <- "Decreased_RNA_Expression" #14

write.csv(RNA_Lip_DESeq_auto_results, "RNA_Lip_DESeq_results_autosomes_combined_sex.csv")

#Autosome Interaction-Term

#Step 3: Construct DESeqDataSet object, removed "instrument" from design because confound perfectly with GSE_number
RNA_Lip_DESeq_auto_IT <- DESeqDataSetFromMatrix(countData = filt_read_counts_autosomes, colData = RNA_Metadata_ex_rem, design = ~ GSE_number + predicted_fetal_sex*disease_group)

#Pre-filtering: removing rows with low gene counts, keeping rows that have at least 10 reads adding across all samples
pre_filter <- rowSums(counts(RNA_Lip_DESeq_auto_IT)) >=10
RNA_Lip_DESeq_auto_IT <- RNA_Lip_DESeq_auto_IT[pre_filter,]

#Set the factor level
RNA_Lip_DESeq_auto_IT$disease_group <- relevel(RNA_Lip_DESeq_auto_IT$disease_group, ref = "Control")

#Step 4: Run DESeq
RNA_Lip_DESeq_auto_IT <- DESeq(RNA_Lip_DESeq_auto_IT)
## Positive LFC means gene went up more in males, negative LFC means gene went up more in females
RNA_Lip_DESeq_auto_IT_results <- results(RNA_Lip_DESeq_auto_IT, name="predicted_fetal_sexM.disease_groupPE")
RNA_Lip_DESeq_auto_IT_results <- as.data.frame(RNA_Lip_DESeq_auto_IT_results)
RNA_Lip_DESeq_auto_IT_results <- merge(RNA_Lip_DESeq_auto_IT_results, lipid_ensembl_list[, c("gene_id", "gene_symbol", "seqnames")], by.x = "row.names", by.y = "gene_id")
RNA_Lip_DESeq_auto_IT_results <- RNA_Lip_DESeq_auto_IT_results[!is.na(RNA_Lip_DESeq_auto_IT_results$padj),]

RNA_Lip_DESeq_auto_IT_results$Expression_Status <- "Not_Biologically_Significant" #4753
RNA_Lip_DESeq_auto_IT_results$Expression_Status[RNA_Lip_DESeq_auto_IT_results$log2FoldChange > 0.00 & RNA_Lip_DESeq_auto_IT_results$padj <0.05] <- "Trending_Towards_Increased_RNA_Expression" #0
RNA_Lip_DESeq_auto_IT_results$Expression_Status[RNA_Lip_DESeq_auto_IT_results$log2FoldChange < 0.00 & RNA_Lip_DESeq_auto_IT_results$padj <0.05] <- "Trending_Towards_Decreased_RNA_Expression" #0
RNA_Lip_DESeq_auto_IT_results$Expression_Status[RNA_Lip_DESeq_auto_IT_results$log2FoldChange > 1.00 & RNA_Lip_DESeq_auto_IT_results$padj <0.05] <- "Increased_RNA_Expression" #7
RNA_Lip_DESeq_auto_IT_results$Expression_Status[RNA_Lip_DESeq_auto_IT_results$log2FoldChange < -1.00 & RNA_Lip_DESeq_auto_IT_results$padj <0.05] <- "Decreased_RNA_Expression" #7

write.csv(RNA_Lip_DESeq_auto_IT_results, "RNA_Lip_DESeq_results_autosomes_combined_sex_IT.csv")

#------------------------------------------------

#Female-Specific Differential Expression Analysis, autosomes
Female_Samples <- RNA_Metadata_ex_rem[RNA_Metadata_ex_rem$predicted_fetal_sex == "F",]
filt_read_counts_F<- filt_read_counts[,colnames(filt_read_counts) %in% rownames(Female_Samples)]

filt_read_counts_auto_F <- filt_read_counts_F[rownames(filt_read_counts_F) %in% autosomal_lipids$gene_id,]
filt_read_counts_chrX_F <- filt_read_counts_F[rownames(filt_read_counts_F) %in% chrX_lipids$gene_id,]

all(colnames(filt_read_counts_auto_F) %in% rownames(Female_Samples)) #TRUE
all(colnames(filt_read_counts_auto_F) == rownames(Female_Samples)) #TRUE

RNA_WG_DESeq_auto_F <- DESeqDataSetFromMatrix(countData = filt_read_counts_auto_F, colData = Female_Samples, design = ~ GSE_number + disease_group)

pre_filter <- rowSums(counts(RNA_WG_DESeq_auto_F)) >=10
RNA_WG_DESeq_auto_F <- RNA_WG_DESeq_auto_F[pre_filter,]
RNA_Lip_DESeq_auto_F <- RNA_WG_DESeq_auto_F[rownames(RNA_WG_DESeq_auto_F) %in% lipid_ensembl_list$gene_id]

RNA_Lip_DESeq_auto_F$disease_group <- relevel(RNA_Lip_DESeq_auto_F$disease_group, ref = "Control")
RNA_Lip_DESeq_auto_F <- DESeq(RNA_Lip_DESeq_auto_F)

RNA_Lip_DESeq_auto_F_results <- results(RNA_Lip_DESeq_auto_F)
RNA_Lip_DESeq_auto_F_results <- as.data.frame(RNA_Lip_DESeq_auto_F_results)
RNA_Lip_DESeq_auto_F_results <- merge(RNA_Lip_DESeq_auto_F_results, lipid_ensembl_list[, c("gene_id", "gene_symbol", "seqnames")], by.x = "row.names", by.y = "gene_id")
RNA_Lip_DESeq_auto_F_results <- RNA_Lip_DESeq_auto_F_results[!is.na(RNA_Lip_DESeq_auto_F_results$padj),]

#Labels 
RNA_Lip_DESeq_auto_F_results$Expression_Status <- "Not_Biologically_Significant" #4896
RNA_Lip_DESeq_auto_F_results$Expression_Status[RNA_Lip_DESeq_auto_F_results$log2FoldChange > 0.00 & RNA_Lip_DESeq_auto_F_results$padj <0.05] <- "Trending_Towards_Increased_RNA_Expression" #0
RNA_Lip_DESeq_auto_F_results$Expression_Status[RNA_Lip_DESeq_auto_F_results$log2FoldChange < 0.00 & RNA_Lip_DESeq_auto_F_results$padj <0.05] <- "Trending_Towards_Decreased_RNA_Expression" #0
RNA_Lip_DESeq_auto_F_results$Expression_Status[RNA_Lip_DESeq_auto_F_results$log2FoldChange > 1.00 & RNA_Lip_DESeq_auto_F_results$padj <0.05] <- "Increased_RNA_Expression" #3
RNA_Lip_DESeq_auto_F_results$Expression_Status[RNA_Lip_DESeq_auto_F_results$log2FoldChange < -1.00 & RNA_Lip_DESeq_auto_F_results$padj <0.05] <- "Decreased_RNA_Expression" #1

write.csv(RNA_Lip_DESeq_auto_F_results, "RNA_Lip_DESeq_results_autosomes_F.csv")


#Female-Specific Differential Expression Analysis, chrX

Female_Samples <- RNA_Metadata_ex_rem[RNA_Metadata_ex_rem$predicted_fetal_sex == "F",]

filt_read_counts_chrX_F <- filt_read_counts_chrX[, (colnames(filt_read_counts_chrX) %in% rownames(Female_Samples))]

all(colnames(filt_read_counts_chrX_F) %in% rownames(Female_Samples)) #TRUE
all(colnames(filt_read_counts_chrX_F) == rownames(Female_Samples)) #TRUE

RNA_WG_DESeq_chrX_F <- DESeqDataSetFromMatrix(countData = filt_read_counts_chrX_F, colData = Female_Samples, design = ~ GSE_number + disease_group)

pre_filter <- rowSums(counts(RNA_WG_DESeq_chrX_F)) >=10
RNA_WG_DESeq_chrX_F <- RNA_WG_DESeq_chrX_F[pre_filter,]
RNA_Lip_DESeq_chrX_F <- RNA_WG_DESeq_chrX_F[rownames(RNA_WG_DESeq_chrX_F) %in% lipid_ensembl_list$gene_id]

RNA_Lip_DESeq_chrX_F$disease_group <- relevel(RNA_Lip_DESeq_chrX_F$disease_group, ref = "Control")
RNA_Lip_DESeq_chrX_F <- DESeq(RNA_Lip_DESeq_chrX_F)

RNA_Lip_DESeq_chrX_F_results <- results(RNA_Lip_DESeq_chrX_F)
RNA_Lip_DESeq_chrX_F_results <- as.data.frame(RNA_Lip_DESeq_chrX_F_results)
RNA_Lip_DESeq_chrX_F_results <- merge(RNA_Lip_DESeq_chrX_F_results, lipid_ensembl_list[, c("gene_id", "gene_symbol", "seqnames")], by.x = "row.names", by.y = "gene_id")
RNA_Lip_DESeq_chrX_F_results <- RNA_Lip_DESeq_chrX_F_results[!is.na(RNA_Lip_DESeq_chrX_F_results$padj),]

RNA_Lip_DESeq_chrX_F_results$Expression_Status <- "Not_Biologically_Significant" #147
RNA_Lip_DESeq_chrX_F_results$Expression_Status[RNA_Lip_DESeq_chrX_F_results$log2FoldChange > 0.00 & RNA_Lip_DESeq_chrX_F_results$padj <0.05] <- "Trending_Towards_Increased_RNA_Expression" #0
RNA_Lip_DESeq_chrX_F_results$Expression_Status[RNA_Lip_DESeq_chrX_F_results$log2FoldChange < 0.00 & RNA_Lip_DESeq_chrX_F_results$padj <0.05] <- "Trending_Towards_Decreased_RNA_Expression" #0
RNA_Lip_DESeq_chrX_F_results$Expression_Status[RNA_Lip_DESeq_chrX_F_results$log2FoldChange > 1.00 & RNA_Lip_DESeq_chrX_F_results$padj <0.05] <- "Increased_RNA_Expression" #0
RNA_Lip_DESeq_chrX_F_results$Expression_Status[RNA_Lip_DESeq_chrX_F_results$log2FoldChange < -1.00 & RNA_Lip_DESeq_chrX_F_results$padj <0.05] <- "Decreased_RNA_Expression" #0

write.csv(RNA_Lip_DESeq_chrX_F_results, "RNA_Lip_DESeq_results_chrX_F.csv")

#------------------------------------------------

#Male-Specific Differential Expression Analysis, autosomes

Male_Samples <- RNA_Metadata_ex_rem[RNA_Metadata_ex_rem$predicted_fetal_sex == "M",]

filt_read_counts_auto_M <- filt_read_counts_autosomes[, (colnames(filt_read_counts_autosomes) %in% rownames(Male_Samples))]

all(colnames(filt_read_counts_auto_M) %in% rownames(Male_Samples)) #TRUE
all(colnames(filt_read_counts_auto_M) == rownames(Male_Samples)) #TRUE

RNA_WG_DESeq_auto_M <- DESeqDataSetFromMatrix(countData = filt_read_counts_auto_M, colData = Male_Samples, design = ~ GSE_number + disease_group)

pre_filter <- rowSums(counts(RNA_WG_DESeq_auto_M)) >=10
RNA_WG_DESeq_auto_M <- RNA_WG_DESeq_auto_M[pre_filter,]
RNA_Lip_DESeq_auto_M <- RNA_WG_DESeq_auto_M[rownames(RNA_WG_DESeq_auto_M) %in% lipid_ensembl_list$gene_id]

RNA_Lip_DESeq_auto_M$disease_group <- relevel(RNA_Lip_DESeq_auto_M$disease_group, ref = "Control")
RNA_Lip_DESeq_auto_M <- DESeq(RNA_Lip_DESeq_auto_M)

RNA_Lip_DESeq_auto_M_results <- results(RNA_Lip_DESeq_auto_M)
RNA_Lip_DESeq_auto_M_results <- as.data.frame(RNA_Lip_DESeq_auto_M_results)
RNA_Lip_DESeq_auto_M_results <- merge(RNA_Lip_DESeq_auto_M_results, lipid_ensembl_list[, c("gene_id", "gene_symbol", "seqnames")], by.x = "row.names", by.y = "gene_id")
RNA_Lip_DESeq_auto_M_results <- RNA_Lip_DESeq_auto_M_results[!is.na(RNA_Lip_DESeq_auto_M_results$padj),]

RNA_Lip_DESeq_auto_M_results$Expression_Status <- "Not_Biologically_Significant" #4811
RNA_Lip_DESeq_auto_M_results$Expression_Status[RNA_Lip_DESeq_auto_M_results$log2FoldChange > 0.00 & RNA_Lip_DESeq_auto_M_results$padj <0.05] <- "Trending_Towards_Increased_RNA_Expression" #0
RNA_Lip_DESeq_auto_M_results$Expression_Status[RNA_Lip_DESeq_auto_M_results$log2FoldChange < 0.00 & RNA_Lip_DESeq_auto_M_results$padj <0.05] <- "Trending_Towards_Decreased_RNA_Expression" #0
RNA_Lip_DESeq_auto_M_results$Expression_Status[RNA_Lip_DESeq_auto_M_results$log2FoldChange > 1.00 & RNA_Lip_DESeq_auto_M_results$padj <0.05] <- "Increased_RNA_Expression" #2
RNA_Lip_DESeq_auto_M_results$Expression_Status[RNA_Lip_DESeq_auto_M_results$log2FoldChange < -1.00 & RNA_Lip_DESeq_auto_M_results$padj <0.05] <- "Decreased_RNA_Expression" #0

write.csv(RNA_Lip_DESeq_auto_M_results, "RNA_Lip_DESeq_results_autosomes_M.csv")


#Male-Specific Differential Expression Analysis, chrX

Male_Samples <- RNA_Metadata_ex_rem[RNA_Metadata_ex_rem$predicted_fetal_sex == "M",]

filt_read_counts_chrX_M <- filt_read_counts_chrX[, (colnames(filt_read_counts_chrX) %in% rownames(Male_Samples))]

all(colnames(filt_read_counts_chrX_M) %in% rownames(Male_Samples)) #TRUE
all(colnames(filt_read_counts_chrX_M) == rownames(Male_Samples)) #TRUE

RNA_WG_DESeq_chrX_M <- DESeqDataSetFromMatrix(countData = filt_read_counts_chrX_M, colData = Male_Samples, design = ~ GSE_number + disease_group)

pre_filter <- rowSums(counts(RNA_WG_DESeq_chrX_M)) >=10
RNA_WG_DESeq_chrX_M <- RNA_WG_DESeq_chrX_M[pre_filter,]
RNA_Lip_DESeq_chrX_M <- RNA_WG_DESeq_chrX_M[rownames(RNA_WG_DESeq_chrX_M) %in% lipid_ensembl_list$gene_id]

RNA_Lip_DESeq_chrX_M$disease_group <- relevel(RNA_Lip_DESeq_chrX_M$disease_group, ref = "Control")
RNA_Lip_DESeq_chrX_M <- DESeq(RNA_Lip_DESeq_chrX_M)

RNA_Lip_DESeq_chrX_M_results <- results(RNA_Lip_DESeq_chrX_M)
RNA_Lip_DESeq_chrX_M_results <- as.data.frame(RNA_Lip_DESeq_chrX_M_results)
RNA_Lip_DESeq_chrX_M_results <- merge(RNA_Lip_DESeq_chrX_M_results, lipid_ensembl_list[, c("gene_id", "gene_symbol", "seqnames")], by.x = "row.names", by.y = "gene_id")
RNA_Lip_DESeq_chrX_M_results <- RNA_Lip_DESeq_chrX_M_results[!is.na(RNA_Lip_DESeq_chrX_M_results$padj),]

RNA_Lip_DESeq_chrX_M_results$Expression_Status <- "Not_Biologically_Significant" #144
RNA_Lip_DESeq_chrX_M_results$Expression_Status[RNA_Lip_DESeq_chrX_M_results$log2FoldChange > 0.00 & RNA_Lip_DESeq_chrX_M_results$padj <0.05] <- "Trending_Towards_Increased_RNA_Expression" #0
RNA_Lip_DESeq_chrX_M_results$Expression_Status[RNA_Lip_DESeq_chrX_M_results$log2FoldChange < 0.00 & RNA_Lip_DESeq_chrX_M_results$padj <0.05] <- "Trending_Towards_Decreased_RNA_Expression" #0
RNA_Lip_DESeq_chrX_M_results$Expression_Status[RNA_Lip_DESeq_chrX_M_results$log2FoldChange > 1.00 & RNA_Lip_DESeq_chrX_M_results$padj <0.05] <- "Increased_RNA_Expression" #0
RNA_Lip_DESeq_chrX_M_results$Expression_Status[RNA_Lip_DESeq_chrX_M_results$log2FoldChange < -1.00 & RNA_Lip_DESeq_chrX_M_results$padj <0.05] <- "Decreased_RNA_Expression" #0

write.csv(RNA_Lip_DESeq_chrX_M_results, "RNA_Lip_DESeq_results_chrX_M.csv")


#Male-Specific Differential Expression Analysis, chrY

Male_Samples <- RNA_Metadata_ex_rem[RNA_Metadata_ex_rem$predicted_fetal_sex == "M",]

filt_read_counts_chrY_M <- filt_read_counts_chrY[, (colnames(filt_read_counts_chrY) %in% rownames(Male_Samples))]

all(colnames(filt_read_counts_chrY_M) %in% rownames(Male_Samples)) #TRUE
all(colnames(filt_read_counts_chrY_M) == rownames(Male_Samples)) #TRUE

RNA_WG_DESeq_chrY_M <- DESeqDataSetFromMatrix(countData = filt_read_counts_chrY_M, colData = Male_Samples, design = ~ GSE_number + disease_group)

pre_filter <- rowSums(counts(RNA_WG_DESeq_chrY_M)) >=10
RNA_WG_DESeq_chrY_M <- RNA_WG_DESeq_chrY_M[pre_filter,]
RNA_Lip_DESeq_chrY_M <- RNA_WG_DESeq_chrY_M[rownames(RNA_WG_DESeq_chrY_M) %in% lipid_ensembl_list$gene_id]

RNA_Lip_DESeq_chrY_M$disease_group <- relevel(RNA_Lip_DESeq_chrY_M$disease_group, ref = "Control")
RNA_Lip_DESeq_chrY_M <- DESeq(RNA_Lip_DESeq_chrY_M)

#all genes have equal values for all samples. will not be able to perform differential analysis

# RNA_Lip_DESeq_chrY_M_results <- results(RNA_Lip_DESeq_chrY_M)
# RNA_Lip_DESeq_chrY_M_results <- as.data.frame(RNA_Lip_DESeq_chrY_M_results)
# RNA_Lip_DESeq_chrY_M_results <- merge(RNA_Lip_DESeq_chrY_M_results, lipid_ensembl_list[, c("gene_id", "gene_symbol", "seqnames")], by.x = "row.names", by.y = "gene_id")
# RNA_Lip_DESeq_chrY_M_results <- RNA_Lip_DESeq_chrY_M_results[!is.na(RNA_Lip_DESeq_chrY_M_results$padj),]

# RNA_Lip_DESeq_chrY_M_results$Expression_Status <- "Not_Biologically_Significant"
# RNA_Lip_DESeq_chrY_M_results$Expression_Status[RNA_Lip_DESeq_chrY_M_results$log2FoldChange > 0.00 & RNA_Lip_DESeq_chrY_M_results$padj <0.05] <- "Trending_Towards_Increased_RNA_Expression"
# RNA_Lip_DESeq_chrY_M_results$Expression_Status[RNA_Lip_DESeq_chrY_M_results$log2FoldChange < 0.00 & RNA_Lip_DESeq_chrY_M_results$padj <0.05] <- "Trending_Towards_Decreased_RNA_Expression"
# RNA_Lip_DESeq_chrY_M_results$Expression_Status[RNA_Lip_DESeq_chrY_M_results$log2FoldChange > 1.00 & RNA_Lip_DESeq_chrY_M_results$padj <0.05] <- "Increased_RNA_Expression"
# RNA_Lip_DESeq_chrY_M_results$Expression_Status[RNA_Lip_DESeq_chrY_M_results$log2FoldChange < -1.00 & RNA_Lip_DESeq_chrY_M_results$padj <0.05] <- "Decreased_RNA_Expression"


#write.csv(RNA_Lip_DESeq_chrY_M_results, "RNA_Lip_DESeq_results_chrY_M.csv")


#------------------------------------------------


#Volcano Plots
library(ggplot2) #version 3.5.1
library(gridExtra) #version 2.3
library(ggrepel) #version 0.9.5

#"grey" (no change in expression), "#d02670"- (pink-Increased Expression), "#8a00c4"- (purple-Decreased Expression)

RNA_Lip_DESeq_auto_results$siglabel <- ifelse(RNA_Lip_DESeq_auto_results$Expression_Status %in% c("Increased_RNA_Expression", "Decreased_RNA_Expression"), RNA_Lip_DESeq_auto_results$gene_symbol, NA)
combined_auto <- ggplot(data = RNA_Lip_DESeq_auto_results, aes(x = log2FoldChange, y = -log10(padj), col = Expression_Status)) + 
  geom_point(shape = 19, alpha = 0.3, size = 3) + 
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 18)) +
  ylab("-log10(FDR)") +
  xlab("log2FoldChange") + 
  scale_y_continuous(breaks = seq(0, 17, by = 1), limits = c(0, 17)) +
  scale_x_continuous(breaks = seq(-6, 6, by = 1), limits = c(-6, 6)) +
  scale_color_manual(values = c("Decreased_RNA_Expression" = "#8a00c4", "Increased_RNA_Expression" = "#d02670", "Trending_Towards_Decreased_RNA_Expression" = "grey", "Trending_Towards_Increased_RNA_Expression" = "grey", "Not_Biologically_Significant" = "grey"),
                     guide = "none") +
  geom_vline(xintercept = c(-1, 1), col = "black", linetype = "dashed", linewidth = 0.75) +
  geom_hline(yintercept = c(-log10(0.05)), col = "black", linetype = "dashed", linewidth = 0.75)
  #geom_text_repel(aes(label=siglabel), na.rm = TRUE, max.overlaps = Inf, size = 4, segment.colour = 'grey50')
  
RNA_Lip_DESeq_auto_F_results$siglabel <- ifelse(RNA_Lip_DESeq_auto_F_results$Expression_Status %in% c("Increased_RNA_Expression", "Decreased_RNA_Expression"), RNA_Lip_DESeq_auto_F_results$gene_symbol, NA)
female_auto <- ggplot(data = RNA_Lip_DESeq_auto_F_results, aes(x = log2FoldChange, y = -log10(padj), col = Expression_Status)) + 
  geom_point(shape = 19, alpha = 0.3, size = 3) + 
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 18)) +
  ylab("") +
  xlab("log2FoldChange") + 
  scale_y_continuous(breaks = seq(0, 17, by = 1), limits = c(0, 17)) +
  scale_x_continuous(breaks = seq(-6, 6, by = 1), limits = c(-6, 6)) +
  scale_color_manual(values = c("Decreased_RNA_Expression" = "#8a00c4", "Increased_RNA_Expression" = "#d02670", "Trending_Towards_Decreased_RNA_Expression" = "grey", "Trending_Towards_Increased_RNA_Expression" = "grey", "Not_Biologically_Significant" = "grey"),
                     guide = "none") +
  geom_vline(xintercept = c(-1, 1), col = "black", linetype = "dashed", linewidth = 0.75) +
  geom_hline(yintercept = c(-log10(0.05)), col = "black", linetype = "dashed", linewidth = 0.75)
  #geom_text_repel(aes(label=siglabel), na.rm = TRUE, max.overlaps = Inf, size = 4, segment.colour = 'grey50')


RNA_Lip_DESeq_auto_M_results$siglabel <- ifelse(RNA_Lip_DESeq_auto_M_results$Expression_Status %in% c("Increased_RNA_Expression", "Decreased_RNA_Expression"), RNA_Lip_DESeq_auto_M_results$gene_symbol, NA)
male_auto <- ggplot(data = RNA_Lip_DESeq_auto_M_results, aes(x = log2FoldChange, y = -log10(padj), col = Expression_Status)) + 
  geom_point(shape = 19, alpha = 0.3, size = 3) + 
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 18)) +
  ylab("") +
  xlab("log2FoldChange") + 
  scale_y_continuous(breaks = seq(0, 17, by = 1), limits = c(0, 17)) +
  scale_x_continuous(breaks = seq(-6, 6, by = 1), limits = c(-6, 6)) +
  scale_color_manual(values = c("Decreased_RNA_Expression" = "#8a00c4", "Increased_RNA_Expression" = "#d02670", "Trending_Towards_Decreased_RNA_Expression" = "grey", "Trending_Towards_Increased_RNA_Expression" = "grey", "Not_Biologically_Significant" = "grey"),
                     guide = "none") +
  geom_vline(xintercept = c(-1, 1), col = "black", linetype = "dashed", linewidth = 0.75) +
  geom_hline(yintercept = c(-log10(0.05)), col = "black", linetype = "dashed", linewidth = 0.75) 
  #geom_text_repel(aes(label=siglabel), na.rm = TRUE, max.overlaps = Inf, size = 4, segment.colour = 'grey50')


RNA_Lip_DESeq_chrX_F_results$siglabel <- ifelse(RNA_Lip_DESeq_chrX_F_results$Expression_Status %in% c("Increased_RNA_Expression", "Decreased_RNA_Expression"), RNA_Lip_DESeq_chrX_F_results$gene_symbol, NA)
female_X <- ggplot(data = RNA_Lip_DESeq_chrX_F_results, aes(x = log2FoldChange, y = -log10(padj), col = Expression_Status)) + 
  geom_point(shape = 19, alpha = 0.3, size = 3) + 
  theme_bw() +
  theme(axis.text = element_text(size = 12.5),
        axis.title = element_text(size = 14)) +
  ylab("-log10(FDR)") +
  xlab("log2FoldChange") + 
  scale_y_continuous(breaks = seq(0, 13, by = 1), limits = c(0, 13)) +
  scale_x_continuous(breaks = seq(-6, 6, by = 1), limits = c(-6, 6)) +
  scale_color_manual(values = c("Decreased_RNA_Expression" = "#8a00c4", "Increased_RNA_Expression" = "#d02670", "Trending_Towards_Decreased_RNA_Expression" = "grey", "Trending_Towards_Increased_RNA_Expression" = "grey", "Not_Biologically_Significant" = "grey"),
                     guide = "none") +
  geom_vline(xintercept = c(-1, 1), col = "black", linetype = "dashed", linewidth = 0.75) +
  geom_hline(yintercept = c(-log10(0.05)), col = "black", linetype = "dashed", linewidth = 0.75) +
  geom_text_repel(aes(label=siglabel), na.rm = TRUE, max.overlaps = Inf, size = 4, segment.colour = 'grey50')


RNA_Lip_DESeq_chrX_M_results$siglabel <- ifelse(RNA_Lip_DESeq_chrX_M_results$Expression_Status %in% c("Increased_RNA_Expression", "Decreased_RNA_Expression"), RNA_Lip_DESeq_chrX_M_results$gene_symbol, NA)
male_X <- ggplot(data = RNA_Lip_DESeq_chrX_M_results, aes(x = log2FoldChange, y = -log10(padj), col = Expression_Status)) + 
  geom_point(shape = 19, alpha = 0.3, size = 3) + 
  theme_bw() +
  theme(axis.text = element_text(size = 12.5),
        axis.title = element_text(size = 14)) +
  ylab("") +
  xlab("log2FoldChange") + 
  scale_y_continuous(breaks = seq(0, 13, by = 1), limits = c(0, 13)) +
  scale_x_continuous(breaks = seq(-6, 6, by = 1), limits = c(-6, 6)) +
  scale_color_manual(values = c("Decreased_RNA_Expression" = "#8a00c4", "Increased_RNA_Expression" = "#d02670", "Trending_Towards_Decreased_RNA_Expression" = "grey", "Trending_Towards_Increased_RNA_Expression" = "grey", "Not_Biologically_Significant" = "grey"),
                     guide = "none") +
  geom_vline(xintercept = c(-1, 1), col = "black", linetype = "dashed", linewidth = 0.75) +
  geom_hline(yintercept = c(-log10(0.05)), col = "black", linetype = "dashed", linewidth = 0.75) +
  geom_text_repel(aes(label=siglabel), na.rm = TRUE, max.overlaps = Inf, size = 4, segment.colour = 'grey50')


# RNA_Lip_DESeq_chrY_M_results$siglabel <- ifelse(RNA_Lip_DESeq_chrY_M_results$Expression_Status %in% c("Increased_RNA_Expression", "Decreased_RNA_Expression"), RNA_Lip_DESeq_chrY_M_results$gene_symbol, NA)
# male_Y <- ggplot(data = RNA_Lip_DESeq_chrY_M_results, aes(x = log2FoldChange, y = -log10(padj), col = Expression_Status)) + 
#   geom_vline(xintercept = c(-0.05,0.05), col = "black", linetype = "dashed", linewidth = 0.75) +
#   geom_hline(yintercept = c(-log10(0.05)), col = "black", linetype = "dashed", linewidth = 0.75) +
#   geom_point(shape = 19, alpha = 0.3, size = 3) + 
#   theme_bw() +
#   theme(axis.text = element_text(size = 12.5),
#         axis.title = element_text(size = 14)) +
#   ylab("-log10(FDR)") +
#   xlab("log2FoldChange") + 
#   scale_y_continuous(breaks = seq(0, 6, by = 1), limits = c(0, 6)) +
#   scale_x_continuous(breaks = seq(-2.0, 2.0, by = 0.1), limits = c(-2.0, 2.0)) +
#   scale_color_manual(values = c("#8a00c4","#d02670","grey","grey", "grey"),
#                      guide = "none")

png("./all_autosome_panel_vol_RNAexp_post.png", height = 6, width = 20, units = "in", res = 300)
grid.arrange(combined_auto, female_auto, male_auto, nrow = 1)
dev.off()


png("./combined_autosome_vol_RNAexp.png", height = 9, width = 15, units = "in", res = 300)
grid.arrange(combined_auto, nrow = 1)
#Warning messages:
#   1: Removed 4 rows containing missing values or values outside the scale range
# (`geom_point()`).
dev.off()

png("./fetalsex_autosome_vol_RNAexp_panel.png", height = 9, width = 15, units = "in", res = 300)
grid.arrange(female_auto, male_auto, nrow = 1)
dev.off()

png("./X_vol_RNAexp_panel.png", height = 9, width = 15, units = "in", res = 300)
grid.arrange(female_X, male_X, nrow = 1)
dev.off()





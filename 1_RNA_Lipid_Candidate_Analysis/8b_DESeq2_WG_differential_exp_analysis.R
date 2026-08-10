
#Step 1: Match Lipid-related gene symbols to ensembl IDs using gencode file
BiocManager::install("rtracklayer")
BiocManager::install("data.table")
BiocManager::install("DESeq2")
library(rtracklayer) #version 1.72.0
library(data.table) #version 1.18.4
library(dplyr) #version 1.2.1
library(tidyverse) #version 2.0.0
library(DESeq2) #version 1.52.0

setwd("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_WG_rerun/")

gtf_file <- "../2025_RNA_Lipid_Candidate/genome_mapping/gencode.v48.primary_assembly.annotation.gtf"
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

#Step 3: Construct DESeqDataSet object
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
RNA_WG_DESeq_auto_results$Expression_Status <- "Not_Biologically_Significant" #25726
RNA_WG_DESeq_auto_results$Expression_Status[RNA_WG_DESeq_auto_results$log2FoldChange > 0.00 & RNA_WG_DESeq_auto_results$padj <0.05] <- "Trending_Towards_Increased_RNA_Expression" #1758
RNA_WG_DESeq_auto_results$Expression_Status[RNA_WG_DESeq_auto_results$log2FoldChange < 0.00 & RNA_WG_DESeq_auto_results$padj <0.05] <- "Trending_Towards_Decreased_RNA_Expression" #1169
RNA_WG_DESeq_auto_results$Expression_Status[RNA_WG_DESeq_auto_results$log2FoldChange > 1.00 & RNA_WG_DESeq_auto_results$padj <0.05] <- "Increased_RNA_Expression" #450
RNA_WG_DESeq_auto_results$Expression_Status[RNA_WG_DESeq_auto_results$log2FoldChange < -1.00 & RNA_WG_DESeq_auto_results$padj <0.05] <- "Decreased_RNA_Expression" #148

write.csv(RNA_WG_DESeq_auto_results, "RNA_WG_DESeq_results_autosomes_combined_sex.csv")

#Autosome Interaction-Term

#Step 3: Construct DESeqDataSet object, removed "instrument" from design because confound perfectly with GSE_number
RNA_WG_DESeq_auto_IT <- DESeqDataSetFromMatrix(countData = filt_read_counts_autosomes, colData = RNA_Metadata_ex_rem, design = ~ GSE_number + predicted_fetal_sex*disease_group)

#Pre-filtering: removing rows with low gene counts, keeping rows that have at least 10 reads adding across all samples
pre_filter <- rowSums(counts(RNA_WG_DESeq_auto_IT)) >=10
RNA_WG_DESeq_auto_IT <- RNA_WG_DESeq_auto_IT[pre_filter,]

#Set the factor level
RNA_WG_DESeq_auto_IT$disease_group <- relevel(RNA_WG_DESeq_auto_IT$disease_group, ref = "Control")

#Step 4: Run DESeq
RNA_WG_DESeq_auto_IT <- DESeq(RNA_WG_DESeq_auto_IT)
## Positive LFC means gene went up more in males, negative LFC means gene went up more in females
RNA_WG_DESeq_auto_IT_results <- results(RNA_WG_DESeq_auto_IT, name="predicted_fetal_sexM.disease_groupPE")
RNA_WG_DESeq_auto_IT_results <- as.data.frame(RNA_WG_DESeq_auto_IT_results)
RNA_WG_DESeq_auto_IT_results <- merge(RNA_WG_DESeq_auto_IT_results, ensembl_list[, c("gene_id", "gene_name", "seqnames")], by.x = "row.names", by.y = "gene_id")
RNA_WG_DESeq_auto_IT_results <- RNA_WG_DESeq_auto_IT_results[!is.na(RNA_WG_DESeq_auto_IT_results$padj),]

RNA_WG_DESeq_auto_IT_results$Expression_Status <- "Not_Biologically_Significant" #4753
RNA_WG_DESeq_auto_IT_results$Expression_Status[RNA_WG_DESeq_auto_IT_results$log2FoldChange > 0.00 & RNA_WG_DESeq_auto_IT_results$padj <0.05] <- "Trending_Towards_Increased_RNA_Expression" #0
RNA_WG_DESeq_auto_IT_results$Expression_Status[RNA_WG_DESeq_auto_IT_results$log2FoldChange < 0.00 & RNA_WG_DESeq_auto_IT_results$padj <0.05] <- "Trending_Towards_Decreased_RNA_Expression" #0
RNA_WG_DESeq_auto_IT_results$Expression_Status[RNA_WG_DESeq_auto_IT_results$log2FoldChange > 1.00 & RNA_WG_DESeq_auto_IT_results$padj <0.05] <- "Increased_RNA_Expression" #7
RNA_WG_DESeq_auto_IT_results$Expression_Status[RNA_WG_DESeq_auto_IT_results$log2FoldChange < -1.00 & RNA_WG_DESeq_auto_IT_results$padj <0.05] <- "Decreased_RNA_Expression" #7

write.csv(RNA_WG_DESeq_auto_IT_results, "RNA_WG_DESeq_results_autosomes_combined_sex_IT.csv")

#------------------------------------------------

#ChrX Interaction-Term

#Step 3: Construct DESeqDataSet object, removed "instrument" from design because confound perfectly with GSE_number
RNA_WG_DESeq_chrX_IT <- DESeqDataSetFromMatrix(countData = filt_read_counts_chrX, colData = RNA_Metadata_ex_rem, design = ~ GSE_number + predicted_fetal_sex*disease_group)

#Pre-filtering: removing rows with low gene counts, keeping rows that have at least 10 reads adding across all samples
pre_filter <- rowSums(counts(RNA_WG_DESeq_chrX_IT)) >=10
RNA_WG_DESeq_chrX_IT <- RNA_WG_DESeq_chrX_IT[pre_filter,]

#Set the factor level
RNA_WG_DESeq_chrX_IT$disease_group <- relevel(RNA_WG_DESeq_chrX_IT$disease_group, ref = "Control")

#Step 4: Run DESeq
RNA_WG_DESeq_chrX_IT <- DESeq(RNA_WG_DESeq_chrX_IT)
## Positive LFC means gene went up more in males, negative LFC means gene went up more in females
RNA_WG_DESeq_chrX_IT_results <- results(RNA_WG_DESeq_chrX_IT, name="predicted_fetal_sexM.disease_groupPE")
RNA_WG_DESeq_chrX_IT_results <- as.data.frame(RNA_WG_DESeq_chrX_IT_results)
RNA_WG_DESeq_chrX_IT_results <- merge(RNA_WG_DESeq_chrX_IT_results, WGid_ensembl_list[, c("gene_id", "gene_symbol", "seqnames")], by.x = "row.names", by.y = "gene_id")
RNA_WG_DESeq_chrX_IT_results <- RNA_WG_DESeq_chrX_IT_results[!is.na(RNA_WG_DESeq_chrX_IT_results$padj),]

RNA_WG_DESeq_chrX_IT_results$Expression_Status <- "Not_Biologically_Significant" #4753
RNA_WG_DESeq_chrX_IT_results$Expression_Status[RNA_WG_DESeq_chrX_IT_results$log2FoldChange > 0.00 & RNA_WG_DESeq_chrX_IT_results$padj <0.05] <- "Trending_Towards_Increased_RNA_Expression" #0
RNA_WG_DESeq_chrX_IT_results$Expression_Status[RNA_WG_DESeq_chrX_IT_results$log2FoldChange < 0.00 & RNA_WG_DESeq_chrX_IT_results$padj <0.05] <- "Trending_Towards_Decreased_RNA_Expression" #0
RNA_WG_DESeq_chrX_IT_results$Expression_Status[RNA_WG_DESeq_chrX_IT_results$log2FoldChange > 1.00 & RNA_WG_DESeq_chrX_IT_results$padj <0.05] <- "Increased_RNA_Expression" #7
RNA_WG_DESeq_chrX_IT_results$Expression_Status[RNA_WG_DESeq_chrX_IT_results$log2FoldChange < -1.00 & RNA_WG_DESeq_chrX_IT_results$padj <0.05] <- "Decreased_RNA_Expression" #7

write.csv(RNA_WG_DESeq_chrX_IT_results, "RNA_WG_DESeq_results_chrX_combined_sex_IT.csv")



#------------------------------------------------

#Female-Specific Differential Expression Analysis, autosomes
Female_Samples <- RNA_Metadata_ex_rem[RNA_Metadata_ex_rem$predicted_fetal_sex == "F",]
rownames(Female_Samples) <- Female_Samples$Run
filt_read_counts_F<- filt_read_counts[,colnames(filt_read_counts) %in% rownames(Female_Samples)]

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
RNA_WG_DESeq_auto_F_results$Expression_Status <- "Not_Biologically_Significant" #27249
RNA_WG_DESeq_auto_F_results$Expression_Status[RNA_WG_DESeq_auto_F_results$log2FoldChange > 0.00 & RNA_WG_DESeq_auto_F_results$padj <0.05] <- "Trending_Towards_Increased_RNA_Expression" #840
RNA_WG_DESeq_auto_F_results$Expression_Status[RNA_WG_DESeq_auto_F_results$log2FoldChange < 0.00 & RNA_WG_DESeq_auto_F_results$padj <0.05] <- "Trending_Towards_Decreased_RNA_Expression" #678
RNA_WG_DESeq_auto_F_results$Expression_Status[RNA_WG_DESeq_auto_F_results$log2FoldChange > 1.00 & RNA_WG_DESeq_auto_F_results$padj <0.05] <- "Increased_RNA_Expression" #336
RNA_WG_DESeq_auto_F_results$Expression_Status[RNA_WG_DESeq_auto_F_results$log2FoldChange < -1.00 & RNA_WG_DESeq_auto_F_results$padj <0.05] <- "Decreased_RNA_Expression" #193

write.csv(RNA_WG_DESeq_auto_F_results, "RNA_WG_DESeq_results_autosomes_F.csv")


#Female-Specific Differential Expression Analysis, chrX

all(colnames(filt_read_counts_chrX_F) %in% rownames(Female_Samples)) #TRUE
all(colnames(filt_read_counts_chrX_F) == rownames(Female_Samples)) #TRUE

RNA_WG_DESeq_chrX_F <- DESeqDataSetFromMatrix(countData = filt_read_counts_chrX_F, colData = Female_Samples, design = ~ GSE_number + disease_group)

pre_filter <- rowSums(counts(RNA_WG_DESeq_chrX_F)) >=10
RNA_WG_DESeq_chrX_F <- RNA_WG_DESeq_chrX_F[pre_filter,]
#RNA_Lip_DESeq_chrX_F <- RNA_WG_DESeq_chrX_F[rownames(RNA_WG_DESeq_chrX_F) %in% lipid_ensembl_list$gene_id]

RNA_WG_DESeq_chrX_F$disease_group <- relevel(RNA_WG_DESeq_chrX_F$disease_group, ref = "Control")
RNA_WG_DESeq_chrX_F <- DESeq(RNA_WG_DESeq_chrX_F)

RNA_WG_DESeq_chrX_F_results <- results(RNA_WG_DESeq_chrX_F)
RNA_WG_DESeq_chrX_F_results <- as.data.frame(RNA_WG_DESeq_chrX_F_results)
RNA_WG_DESeq_chrX_F_results <- merge(RNA_WG_DESeq_chrX_F_results, ensembl_list[, c("gene_id", "gene_name", "seqnames")], by.x = "row.names", by.y = "gene_id")
RNA_WG_DESeq_chrX_F_results <- RNA_WG_DESeq_chrX_F_results[!is.na(RNA_WG_DESeq_chrX_F_results$padj),]

RNA_WG_DESeq_chrX_F_results$Expression_Status <- "Not_Biologically_Significant" #799
RNA_WG_DESeq_chrX_F_results$Expression_Status[RNA_WG_DESeq_chrX_F_results$log2FoldChange > 0.00 & RNA_WG_DESeq_chrX_F_results$padj <0.05] <- "Trending_Towards_Increased_RNA_Expression" #35
RNA_WG_DESeq_chrX_F_results$Expression_Status[RNA_WG_DESeq_chrX_F_results$log2FoldChange < 0.00 & RNA_WG_DESeq_chrX_F_results$padj <0.05] <- "Trending_Towards_Decreased_RNA_Expression" #26
RNA_WG_DESeq_chrX_F_results$Expression_Status[RNA_WG_DESeq_chrX_F_results$log2FoldChange > 1.00 & RNA_WG_DESeq_chrX_F_results$padj <0.05] <- "Increased_RNA_Expression" #12
RNA_WG_DESeq_chrX_F_results$Expression_Status[RNA_WG_DESeq_chrX_F_results$log2FoldChange < -1.00 & RNA_WG_DESeq_chrX_F_results$padj <0.05] <- "Decreased_RNA_Expression" #12

write.csv(RNA_WG_DESeq_chrX_F_results, "RNA_WG_DESeq_results_chrX_F.csv")

#------------------------------------------------

#Male-Specific Differential Expression Analysis, autosomes

Male_Samples <- RNA_Metadata_ex_rem[RNA_Metadata_ex_rem$predicted_fetal_sex == "M",]
rownames(Male_Samples) <- Male_Samples$Run
filt_read_counts_M<- filt_read_counts[,colnames(filt_read_counts) %in% rownames(Male_Samples)]

filt_read_counts_auto_M <- filt_read_counts_M[rownames(filt_read_counts_M) %in% autosomes$gene_id,]
filt_read_counts_chrX_M <- filt_read_counts_M[rownames(filt_read_counts_M) %in% chrX$gene_id,]
filt_read_counts_chrY_M <- filt_read_counts_M[rownames(filt_read_counts_M) %in% chrY$gene_id,]

all(colnames(filt_read_counts_auto_M) %in% rownames(Male_Samples)) #TRUE
all(colnames(filt_read_counts_auto_M) == rownames(Male_Samples)) #TRUE

RNA_WG_DESeq_auto_M <- DESeqDataSetFromMatrix(countData = filt_read_counts_auto_M, colData = Male_Samples, design = ~ GSE_number + disease_group)

pre_filter <- rowSums(counts(RNA_WG_DESeq_auto_M)) >=10
RNA_WG_DESeq_auto_M <- RNA_WG_DESeq_auto_M[pre_filter,]
#RNA_Lip_DESeq_auto_M <- RNA_WG_DESeq_auto_M[rownames(RNA_WG_DESeq_auto_M) %in% lipid_ensembl_list$gene_id]

RNA_WG_DESeq_auto_M$disease_group <- relevel(RNA_WG_DESeq_auto_M$disease_group, ref = "Control")
RNA_WG_DESeq_auto_M <- DESeq(RNA_WG_DESeq_auto_M)

RNA_WG_DESeq_auto_M_results <- results(RNA_WG_DESeq_auto_M)
RNA_WG_DESeq_auto_M_results <- as.data.frame(RNA_WG_DESeq_auto_M_results)
RNA_WG_DESeq_auto_M_results <- merge(RNA_WG_DESeq_auto_M_results, ensembl_list[, c("gene_id", "gene_name", "seqnames")], by.x = "row.names", by.y = "gene_id")
RNA_WG_DESeq_auto_M_results <- RNA_WG_DESeq_auto_M_results[!is.na(RNA_WG_DESeq_auto_M_results$padj),]

RNA_WG_DESeq_auto_M_results$Expression_Status <- "Not_Biologically_Significant" #19419
RNA_WG_DESeq_auto_M_results$Expression_Status[RNA_WG_DESeq_auto_M_results$log2FoldChange > 0.00 & RNA_WG_DESeq_auto_M_results$padj <0.05] <- "Trending_Towards_Increased_RNA_Expression" #212
RNA_WG_DESeq_auto_M_results$Expression_Status[RNA_WG_DESeq_auto_M_results$log2FoldChange < 0.00 & RNA_WG_DESeq_auto_M_results$padj <0.05] <- "Trending_Towards_Decreased_RNA_Expression" #177
RNA_WG_DESeq_auto_M_results$Expression_Status[RNA_WG_DESeq_auto_M_results$log2FoldChange > 1.00 & RNA_WG_DESeq_auto_M_results$padj <0.05] <- "Increased_RNA_Expression" #242
RNA_WG_DESeq_auto_M_results$Expression_Status[RNA_WG_DESeq_auto_M_results$log2FoldChange < -1.00 & RNA_WG_DESeq_auto_M_results$padj <0.05] <- "Decreased_RNA_Expression" #64

write.csv(RNA_WG_DESeq_auto_M_results, "RNA_WG_DESeq_results_autosomes_M.csv")


#Male-Specific Differential Expression Analysis, chrX

all(colnames(filt_read_counts_chrX_M) %in% rownames(Male_Samples)) #TRUE
all(colnames(filt_read_counts_chrX_M) == rownames(Male_Samples)) #TRUE

RNA_WG_DESeq_chrX_M <- DESeqDataSetFromMatrix(countData = filt_read_counts_chrX_M, colData = Male_Samples, design = ~ GSE_number + disease_group)

pre_filter <- rowSums(counts(RNA_WG_DESeq_chrX_M)) >=10
RNA_WG_DESeq_chrX_M <- RNA_WG_DESeq_chrX_M[pre_filter,]
#RNA_Lip_DESeq_chrX_M <- RNA_WG_DESeq_chrX_M[rownames(RNA_WG_DESeq_chrX_M) %in% lipid_ensembl_list$gene_id]

RNA_WG_DESeq_chrX_M$disease_group <- relevel(RNA_WG_DESeq_chrX_M$disease_group, ref = "Control")
RNA_WG_DESeq_chrX_M <- DESeq(RNA_WG_DESeq_chrX_M)

RNA_WG_DESeq_chrX_M_results <- results(RNA_WG_DESeq_chrX_M)
RNA_WG_DESeq_chrX_M_results <- as.data.frame(RNA_WG_DESeq_chrX_M_results)
RNA_WG_DESeq_chrX_M_results <- merge(RNA_WG_DESeq_chrX_M_results, ensembl_list[, c("gene_id", "gene_name", "seqnames")], by.x = "row.names", by.y = "gene_id")
RNA_WG_DESeq_chrX_M_results <- RNA_WG_DESeq_chrX_M_results[!is.na(RNA_WG_DESeq_chrX_M_results$padj),]

RNA_WG_DESeq_chrX_M_results$Expression_Status <- "Not_Biologically_Significant" #765
RNA_WG_DESeq_chrX_M_results$Expression_Status[RNA_WG_DESeq_chrX_M_results$log2FoldChange > 0.00 & RNA_WG_DESeq_chrX_M_results$padj <0.05] <- "Trending_Towards_Increased_RNA_Expression" #7
RNA_WG_DESeq_chrX_M_results$Expression_Status[RNA_WG_DESeq_chrX_M_results$log2FoldChange < 0.00 & RNA_WG_DESeq_chrX_M_results$padj <0.05] <- "Trending_Towards_Decreased_RNA_Expression" #6
RNA_WG_DESeq_chrX_M_results$Expression_Status[RNA_WG_DESeq_chrX_M_results$log2FoldChange > 1.00 & RNA_WG_DESeq_chrX_M_results$padj <0.05] <- "Increased_RNA_Expression" #7
RNA_WG_DESeq_chrX_M_results$Expression_Status[RNA_WG_DESeq_chrX_M_results$log2FoldChange < -1.00 & RNA_WG_DESeq_chrX_M_results$padj <0.05] <- "Decreased_RNA_Expression" #2

write.csv(RNA_WG_DESeq_chrX_M_results, "RNA_WG_DESeq_results_chrX_M.csv")


#Male-Specific Differential Expression Analysis, chrY

Male_Samples <- RNA_Metadata_ex_rem[RNA_Metadata_ex_rem$predicted_fetal_sex == "M",]

filt_read_counts_chrY_M <- filt_read_counts_chrY[, (colnames(filt_read_counts_chrY) %in% rownames(Male_Samples))]

all(colnames(filt_read_counts_chrY_M) %in% rownames(Male_Samples)) #TRUE
all(colnames(filt_read_counts_chrY_M) == rownames(Male_Samples)) #TRUE

RNA_WG_DESeq_chrY_M <- DESeqDataSetFromMatrix(countData = filt_read_counts_chrY_M, colData = Male_Samples, design = ~ GSE_number + disease_group)

pre_filter <- rowSums(counts(RNA_WG_DESeq_chrY_M)) >=10
RNA_WG_DESeq_chrY_M <- RNA_WG_DESeq_chrY_M[pre_filter,]
#RNA_Lip_DESeq_chrY_M <- RNA_WG_DESeq_chrY_M[rownames(RNA_WG_DESeq_chrY_M) %in% lipid_ensembl_list$gene_id]

RNA_WG_DESeq_chrY_M$disease_group <- relevel(RNA_WG_DESeq_chrY_M$disease_group, ref = "Control")
RNA_WG_DESeq_chrY_M <- DESeq(RNA_WG_DESeq_chrY_M)

RNA_WG_DESeq_chrY_M_results <- results(RNA_WG_DESeq_chrY_M)
RNA_WG_DESeq_chrY_M_results <- as.data.frame(RNA_WG_DESeq_chrY_M_results)
RNA_WG_DESeq_chrY_M_results <- merge(RNA_WG_DESeq_chrY_M_results, ensembl_list[, c("gene_id", "gene_name", "seqnames")], by.x = "row.names", by.y = "gene_id")
RNA_WG_DESeq_chrY_M_results <- RNA_WG_DESeq_chrY_M_results[!is.na(RNA_WG_DESeq_chrY_M_results$padj),]

RNA_WG_DESeq_chrY_M_results$Expression_Status <- "Not_Biologically_Significant" #74
RNA_WG_DESeq_chrY_M_results$Expression_Status[RNA_WG_DESeq_chrY_M_results$log2FoldChange > 0.00 & RNA_WG_DESeq_chrY_M_results$padj <0.05] <- "Trending_Towards_Increased_RNA_Expression" #0
RNA_WG_DESeq_chrY_M_results$Expression_Status[RNA_WG_DESeq_chrY_M_results$log2FoldChange < 0.00 & RNA_WG_DESeq_chrY_M_results$padj <0.05] <- "Trending_Towards_Decreased_RNA_Expression" #0
RNA_WG_DESeq_chrY_M_results$Expression_Status[RNA_WG_DESeq_chrY_M_results$log2FoldChange > 1.00 & RNA_WG_DESeq_chrY_M_results$padj <0.05] <- "Increased_RNA_Expression" #0
RNA_WG_DESeq_chrY_M_results$Expression_Status[RNA_WG_DESeq_chrY_M_results$log2FoldChange < -1.00 & RNA_WG_DESeq_chrY_M_results$padj <0.05] <- "Decreased_RNA_Expression" #0

write.csv(RNA_WG_DESeq_chrY_M_results, "RNA_WG_DESeq_results_chrY_M.csv")

#------------------------------------------------

#Volcano Plots
library(ggplot2) #version 3.5.1
library(gridExtra) #version 2.3
library(ggrepel) #version 0.9.5

#"grey" (no change in expression), "#d02670"- (pink-Increased Expression), "#8a00c4"- (purple-Decreased Expression)

RNA_WG_DESeq_auto_results$siglabel <- ifelse(RNA_WG_DESeq_auto_results$Expression_Status %in% c("Increased_RNA_Expression", "Decreased_RNA_Expression"), RNA_WG_DESeq_auto_results$gene_name, NA)
combined_auto <- ggplot(data = RNA_WG_DESeq_auto_results, aes(x = log2FoldChange, y = -log10(padj), col = Expression_Status)) + 
  geom_point(shape = 19, alpha = 0.3, size = 5) + 
  theme_bw() +
  labs(title = "Combined-Sex Autosomal\nWhole-Genome DEGs") +
  theme(plot.title = element_text(size = 17, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18)) +
  ylab("-log10(FDR)") +
  xlab("log2FoldChange") + 
  scale_y_continuous(breaks = seq(0, 19, by = 1), limits = c(0, 19)) +
  scale_x_continuous(breaks = seq(-8, 8, by = 1), limits = c(-8, 8)) +
  scale_color_manual(values = c("#8a00c4","#d02670","grey","grey", "grey"),
                     guide = "none") +
  geom_vline(xintercept = c(-1, 1), col = "black", linetype = "dashed", linewidth = 0.75) +
  geom_hline(yintercept = c(-log10(0.05)), col = "black", linetype = "dashed", linewidth = 0.75)
  #geom_text_repel(aes(label=siglabel), na.rm = TRUE, max.overlaps = Inf, size = 4, segment.colour = 'grey50')
  
RNA_WG_DESeq_auto_F_results$siglabel <- ifelse(RNA_WG_DESeq_auto_F_results$Expression_Status %in% c("Increased_RNA_Expression", "Decreased_RNA_Expression"), RNA_WG_DESeq_auto_F_results$gene_name, NA)
female_auto <- ggplot(data = RNA_WG_DESeq_auto_F_results, aes(x = log2FoldChange, y = -log10(padj), col = Expression_Status)) + 
  geom_point(shape = 19, alpha = 0.3, size = 5) + 
  theme_bw() +
  labs(title = "Fetal Female Autosomal\nWhole-Genome DEGs") +
  theme(plot.title = element_text(size = 17, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18)) +
  ylab("") +
  xlab("log2FoldChange") + 
  scale_y_continuous(breaks = seq(0, 19, by = 1), limits = c(0, 19)) +
  scale_x_continuous(breaks = seq(-8, 8, by = 1), limits = c(-8, 8)) +
  scale_color_manual(values = c("#8a00c4","#d02670","grey","grey", "grey"),
                     guide = "none") +
  geom_vline(xintercept = c(-1, 1), col = "black", linetype = "dashed", linewidth = 0.75) +
  geom_hline(yintercept = c(-log10(0.05)), col = "black", linetype = "dashed", linewidth = 0.75
) 
  #geom_text_repel(aes(label=siglabel), na.rm = TRUE, max.overlaps = Inf, size = 4, segment.colour = 'grey50')


RNA_WG_DESeq_auto_M_results$siglabel <- ifelse(RNA_WG_DESeq_auto_M_results$Expression_Status %in% c("Increased_RNA_Expression", "Decreased_RNA_Expression"), RNA_WG_DESeq_auto_M_results$gene_name, NA)
male_auto <- ggplot(data = RNA_WG_DESeq_auto_M_results, aes(x = log2FoldChange, y = -log10(padj), col = Expression_Status)) + 
  geom_point(shape = 19, alpha = 0.3, size = 5) + 
  theme_bw() +
  labs(title = "Fetal Male Autosomal\nWhole-Genome DEGs") +
  theme(plot.title = element_text(size = 17, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18)) +
  ylab("") +
  xlab("log2FoldChange") + 
  scale_y_continuous(breaks = seq(0, 19, by = 1), limits = c(0, 19)) +
  scale_x_continuous(breaks = seq(-8, 8, by = 1), limits = c(-8, 8)) +
  scale_color_manual(values = c("#8a00c4","#d02670","grey","grey", "grey"),
                     guide = "none") +
  geom_vline(xintercept = c(-1, 1), col = "black", linetype = "dashed", linewidth = 0.75) +
  geom_hline(yintercept = c(-log10(0.05)), col = "black", linetype = "dashed", linewidth = 0.75) 
  #geom_text_repel(aes(label=siglabel), na.rm = TRUE, max.overlaps = Inf, size = 4, segment.colour = 'grey50')


RNA_WG_DESeq_chrX_F_results$siglabel <- ifelse(RNA_WG_DESeq_chrX_F_results$Expression_Status %in% c("Increased_RNA_Expression", "Decreased_RNA_Expression"), RNA_WG_DESeq_chrX_F_results$gene_name, NA)
female_X <- ggplot(data = RNA_WG_DESeq_chrX_F_results, aes(x = log2FoldChange, y = -log10(padj), col = Expression_Status)) + 
  geom_point(shape = 19, alpha = 0.5, size = 5) + 
  theme_bw() +
  labs(title = "Fetal Female X-Chr\nWhole-Genome DEGs") +
  theme(plot.title = element_text(size = 17, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18)) +
  ylab("-log10(FDR)") +
  xlab("log2FoldChange") + 
  scale_y_continuous(breaks = seq(0, 7, by = 1), limits = c(0, 7)) +
  scale_x_continuous(breaks = seq(-6, 6, by = 1), limits = c(-6, 6)) +
  scale_color_manual(values = c("#8a00c4","#d02670","grey","grey", "grey"),
                     guide = "none") +
  geom_vline(xintercept = c(-1, 1), col = "black", linetype = "dashed", linewidth = 0.75) +
  geom_hline(yintercept = c(-log10(0.05)), col = "black", linetype = "dashed", linewidth = 0.75) 
  #geom_text_repel(aes(label=siglabel), na.rm = TRUE, max.overlaps = Inf, size = 4, segment.colour = 'grey50')


RNA_WG_DESeq_chrX_M_results$siglabel <- ifelse(RNA_WG_DESeq_chrX_M_results$Expression_Status %in% c("Increased_RNA_Expression", "Decreased_RNA_Expression"), RNA_WG_DESeq_chrX_M_results$gene_name, NA)
male_X <- ggplot(data = RNA_WG_DESeq_chrX_M_results, aes(x = log2FoldChange, y = -log10(padj), col = Expression_Status)) + 
  geom_point(shape = 19, alpha = 0.5, size = 5) + 
  theme_bw() +
  labs(title = "Fetal Male X-Chr\nWhole-Genome DEGs") +
  theme(plot.title = element_text(size = 17, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18)) +
  ylab("") +
  xlab("log2FoldChange") + 
  scale_y_continuous(breaks = seq(0, 7, by = 1), limits = c(0, 7)) +
  scale_x_continuous(breaks = seq(-6, 6, by = 1), limits = c(-6, 6)) +
  scale_color_manual(values = c("#8a00c4","#d02670","grey","grey", "grey"),
                     guide = "none") +
  geom_vline(xintercept = c(-1, 1), col = "black", linetype = "dashed", linewidth = 0.75) +
  geom_hline(yintercept = c(-log10(0.05)), col = "black", linetype = "dashed", linewidth = 0.75) 
  #geom_text_repel(aes(label=siglabel), na.rm = TRUE, max.overlaps = Inf, size = 4, segment.colour = 'grey50')


# RNA_WG_DESeq_chrY_M_results$siglabel <- ifelse(RNA_WG_DESeq_chrY_M_results$Expression_Status %in% c("Increased_RNA_Expression", "Decreased_RNA_Expression"), RNA_WG_DESeq_chrY_M_results$gene_name, NA)
# male_Y <- ggplot(data = RNA_WG_DESeq_chrY_M_results, aes(x = log2FoldChange, y = -log10(padj), col = Expression_Status)) +
#   geom_point(shape = 19, alpha = 0.3, size = 3) + 
#   theme_bw() +
#   theme(axis.text = element_text(size = 14),
#         axis.title = element_text(size = 18)) +
#   ylab("") +
#   xlab("log2FoldChange") + 
#   scale_y_continuous(breaks = seq(0, 22, by = 1), limits = c(0, 22)) +
#   scale_x_continuous(breaks = seq(-6, 6, by = 1), limits = c(-6, 6)) +
#   geom_vline(xintercept = c(-1, 1), col = "black", linetype = "dashed", linewidth = 0.75) +
#   geom_hline(yintercept = c(-log10(0.05)), col = "black", linetype = "dashed", linewidth = 0.75) +
#   scale_color_manual(values = c("#8a00c4","#d02670","grey","grey", "grey"),
#                      guide = "none")  
#   #geom_text_repel(aes(label=siglabel), na.rm = TRUE, max.overlaps = Inf, size = 4, segment.colour = 'grey50')


png("./WG_panel_autosomes_vol_RNAexp.png", height = 9, width = 15, units = "in", res = 300)
grid.arrange(combined_auto, female_auto, male_auto, nrow = 1)
dev.off()

png("./WG_combined_autosome_vol_RNAexp.png", height = 9, width = 15, units = "in", res = 300)
grid.arrange(combined_auto, nrow = 1)
dev.off()

png("./WG_fetalsex_autosome_vol_RNAexp_panel.png", height = 9, width = 15, units = "in", res = 300)
grid.arrange(female_auto, male_auto, nrow = 1)
dev.off()

png("./WG_chrX_vol_RNAexp_panel.png", height = 9, width = 15, units = "in", res = 300)
grid.arrange(female_X, male_X, nrow = 1)
dev.off()

png("./WG_pval_dist_deseq.png", height = 9, width = 20, units = "in", res = 300)
hist(RNA_WG_DESeq_auto_results$pval, breaks=50, col="skyblue", 
     main="P-value Distribution", xlab="Unadjusted P-value")
dev.off()


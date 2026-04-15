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

cd /workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/DESeq_rerun

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

#Labels 
RNA_Lip_DESeq_auto_results$Expression_Status <- "Not_Biologically_Significant" #4753
RNA_Lip_DESeq_auto_results$Expression_Status[RNA_Lip_DESeq_auto_results$log2FoldChange > 0.00 & RNA_Lip_DESeq_auto_results$padj <0.05] <- "Trending_Towards_Increased_RNA_Expression" #0
RNA_Lip_DESeq_auto_results$Expression_Status[RNA_Lip_DESeq_auto_results$log2FoldChange < 0.00 & RNA_Lip_DESeq_auto_results$padj <0.05] <- "Trending_Towards_Decreased_RNA_Expression" #0
RNA_Lip_DESeq_auto_results$Expression_Status[RNA_Lip_DESeq_auto_results$log2FoldChange > 1.00 & RNA_Lip_DESeq_auto_results$padj <0.05] <- "Increased_RNA_Expression" #7
RNA_Lip_DESeq_auto_results$Expression_Status[RNA_Lip_DESeq_auto_results$log2FoldChange < -1.00 & RNA_Lip_DESeq_auto_results$padj <0.05] <- "Decreased_RNA_Expression" #7

write.csv(RNA_Lip_DESeq_auto_results, "RNA_Lip_DESeq_results_autosomes_combined_sex_IT.csv")


#Identify which differentially methylated genes are common and unique to male and female

##Lipid-Candidate
RNA_Lip_auto_F <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/DESeq_rerun/RNA_Lip_DESeq_results_autosomes_F.csv")
RNA_Lip_auto_M <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/DESeq_rerun/RNA_Lip_DESeq_results_autosomes_M.csv")

affy_Lip_auto_F <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/affymetrix_validation_rerun/affy_results_auto_F.csv")
affy_Lip_auto_M <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/affymetrix_validation_rerun/affy_results_auto_M.csv")


RNA_Lip_auto_F <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/DESeq_rerun/RNA_Lip_DESeq_results_autosomes_F.csv")
RNA_Lip_auto_M <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/DESeq_rerun/RNA_Lip_DESeq_results_autosomes_M.csv")
bio_sig_F <- RNA_Lip_auto_F[RNA_Lip_auto_F$Expression_Status %in% (c("Decreased_RNA_Expression", "Increased_RNA_Expression")), ]
bio_sig_M <- RNA_Lip_auto_M[RNA_Lip_auto_M$Expression_Status %in% (c("Decreased_RNA_Expression", "Increased_RNA_Expression")), ]

#Genes sig only in females and only iin males
only_bio_sig_F <- bio_sig_F[!(bio_sig_F$gene_symbol %in% bio_sig_M$gene_symbol),] #173 genes
only_bio_sig_M <- bio_sig_M[!(bio_sig_M$gene_symbol %in% bio_sig_F$gene_symbol),] #242 genes
#genes in both but with pvals + logFC in females, then males
female_male_DEG_F <- bio_sig_F[(bio_sig_F$gene_symbol %in% bio_sig_M$gene_symbol),] #259 genes
female_male_DEG_M <- bio_sig_M[(bio_sig_M$gene_symbol %in% bio_sig_F$gene_symbol),]

RNA_Lip_auto_F$unique <- "Neither"
RNA_Lip_auto_F$unique[RNA_Lip_auto_F$gene_symbol %in% only_bio_sig_F$gene_symbol] <- "Only_Female"
RNA_Lip_auto_F$unique[RNA_Lip_auto_F$gene_symbol %in% female_male_DEG_F$gene_symbol] <- "Both"

RNA_Lip_auto_M$unique <- "Neither"
RNA_Lip_auto_M$unique[RNA_Lip_auto_M$gene_symbol %in% only_bio_sig_M$gene_symbol] <- "Only_Male"
RNA_Lip_auto_M$unique[RNA_Lip_auto_M$gene_symbol %in% female_male_DEG_M$gene_symbol] <- "Both"

write.csv(RNA_Lip_auto_F, "RNA_Lip_DESeq_results_autosomes_F.csv")
write.csv(RNA_Lip_auto_M, "RNA_Lip_DESeq_results_autosomes_M.csv")

##Whole Genome

RNA_Lip_auto_F_WG <-read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_WG_rerun/RNA_WG_DESeq_results_autosomes_F.csv")
RNA_Lip_auto_M_WG <-read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_WG_rerun/RNA_WG_DESeq_results_autosomes_M.csv")
bio_sig_F <- RNA_Lip_auto_F_WG[RNA_Lip_auto_F_WG$Expression_Status %in% (c("Decreased_RNA_Expression", "Increased_RNA_Expression")), ]
bio_sig_M <- RNA_Lip_auto_M_WG[RNA_Lip_auto_M_WG$Expression_Status %in% (c("Decreased_RNA_Expression", "Increased_RNA_Expression")), ]

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










#Group sig differentially expressed genes into groups to parse out pathways involved in PE sex-differences
#NVM doesn't work because you need the non-sig genes to establish a baseline for ermineJ
RNA_Lip_auto_F <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/DESeq_rerun/RNA_Lip_DESeq_results_autosomes_F.csv")
RNA_Lip_auto_M <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/DESeq_rerun/RNA_Lip_DESeq_results_autosomes_M.csv")
bio_sig_F <- RNA_Lip_auto_F[RNA_Lip_auto_F$Expression_Status %in% (c("Decreased_RNA_Expression", "Increased_RNA_Expression")), ]
bio_sig_M <- RNA_Lip_auto_M[RNA_Lip_auto_M$Expression_Status %in% (c("Decreased_RNA_Expression", "Increased_RNA_Expression")), ]

#Genes sig only in females and only iin males
only_bio_sig_F <- bio_sig_F[!(bio_sig_F$gene_symbol %in% bio_sig_M$gene_symbol),] #173 genes
only_bio_sig_M <- bio_sig_M[!(bio_sig_M$gene_symbol %in% bio_sig_F$gene_symbol),] #242 genes
#genes in both but with pvals + logFC in females, then males
female_male_DEG_F <- bio_sig_F[(bio_sig_F$gene_symbol %in% bio_sig_M$gene_symbol),] #259 genes
female_male_DEG_M <- bio_sig_M[(bio_sig_M$gene_symbol %in% bio_sig_F$gene_symbol),]

#Run Pathway Enrichment on each of 4 groups


 if(!file.exists('Generic_human_ncbiIds_noParents.an.txt.gz')){
    system('wget https://gemma.msl.ubc.ca/annots/Generic_human_ncbiIds_noParents.an.txt.gz')}
NCBI <- read.table('Generic_human_ncbiIds_noParents.an.txt.gz', sep = '\t', header = T, quote = "")

# Female Only
only_bio_sig_F_pval <- only_bio_sig_F[, c("gene_symbol", "pvalue")] # 4710 genes
only_bio_sig_F_pval <- only_bio_sig_F_pval[!is.na(only_bio_sig_F_pval$pvalue),] #4710
only_bio_sig_F_pval <- only_bio_sig_F_pval[!duplicated(only_bio_sig_F_pval$gene_symbol),] # 4708 lipid genes investigated

only_bio_sig_F_reference <- NCBI[NCBI$GeneSymbols %in% only_bio_sig_F_pval$gene_symbol, ] #4698 reference genes

# Male Only
only_bio_sig_M_pval <- only_bio_sig_M[, c("gene_symbol", "pvalue")] # 4614 probes
only_bio_sig_M_pval <- only_bio_sig_M_pval[!is.na(only_bio_sig_M_pval$pvalue),] #4858
only_bio_sig_M_pval <- only_bio_sig_M_pval[!duplicated(only_bio_sig_M_pval$gene_symbol),] # 4613 lipid genes investigated

only_bio_sig_M_reference <- NCBI[NCBI$GeneSymbols %in% only_bio_sig_M_pval$gene_symbol, ] #4607 genes

# Both F pvals Chromosomes
female_male_DEG_F_pval <- female_male_DEG_F[, c("gene_symbol", "pvalue")] # 149 genes
female_male_DEG_F_pval <- female_male_DEG_F_pval[!is.na(female_male_DEG_F_pval$pvalue),] #149
female_male_DEG_F_pval <- female_male_DEG_F_pval[!duplicated(female_male_DEG_F_pval$gene_symbol),] # 149 lipid genes investigated

female_male_DEG_F_reference <- NCBI[NCBI$GeneSymbols %in% female_male_DEG_F_pval$gene_symbol, ] # 149 reference genes

# Both M pvals Chromosome
female_male_DEG_M_pval <- female_male_DEG_M[, c("gene_symbol", "pvalue")] # 147 probes
female_male_DEG_M_pval <- female_male_DEG_M_pval[!is.na(female_male_DEG_M_pval$pvalue),] #147
female_male_DEG_M_pval <- female_male_DEG_M_pval[!duplicated(female_male_DEG_M_pval$gene_symbol),] # 147 lipid genes investigated

female_male_DEG_M_reference <- NCBI[NCBI$GeneSymbols %in% female_male_DEG_M_pval$gene_symbol, ] #147 genes

# Part 2: Pathway Enrichment Analysis

library(ermineR) #version 1.0.3.9000
library(dplyr) #version 1.1.4
library(ggplot2) #version 3.5.1

Sys.setenv('JAVA_HOME' = '/usr/lib/jvm/java-21-openjdk/')

only_F_auto_pr_out <- precRecall(annotation = only_bio_sig_F_reference, 
                    scores = only_bio_sig_F_pval,
                    scoreColumn = 2,
                    logTrans = T,
                    bigIsBetter = F,
                    iterations = 10000)

only_F_auto_pathway_analysis <- as.data.frame(only_F_auto_pr_out$results)

sig_pathways_F_auto <- only_F_auto_pathway_analysis[only_F_auto_pathway_analysis$CorrectedPvalue < 0.05 & only_F_auto_pathway_analysis$Multifunctionality < 0.5,] #0
sig_pathways_F_auto <- sig_pathways_F_auto %>% distinct(GeneMembers, .keep_all = T)
write.csv(sig_pathways_F_auto, '2025_F_auto_enriched_lipid_pathways.csv')

only_M_auto_pr_out <- precRecall(annotation = only_bio_sig_M_reference, 
                    scores = only_bio_sig_M_pval,
                    scoreColumn = 2,
                    logTrans = T,
                    bigIsBetter = F,
                    iterations = 10000)

only_M_auto_pathway_analysis <- as.data.frame(only_M_auto_pr_out$results)

sig_pathways_M_auto <- only_M_auto_pathway_analysis[only_M_auto_pathway_analysis$CorrectedPvalue < 0.05 & only_M_auto_pathway_analysis$Multifunctionality < 0.5,] #0
sig_pathways_M_auto <- sig_pathways_M_auto %>% distinct(GeneMembers, .keep_all = T)
write.csv(sig_pathways_M_auto, '2025_M_auto_enriched_lipid_pathways.csv')

F_M_DEG_F_pr_out <- precRecall(annotation = female_male_DEG_F_reference, 
                    scores = female_male_DEG_F_pval,
                    scoreColumn = 2, 
                    logTrans = T,
                    bigIsBetter = F,
                    iterations = 10000)

F_M_DEG_F_pathway_analysis <- as.data.frame(F_M_DEG_F_pr_out$results)

sig_pathways_F_M_DEG_F <- F_M_DEG_F_pathway_analysis[F_M_DEG_F_pathway_analysis$CorrectedPvalue < 0.05 & F_M_DEG_F_pathway_analysis$Multifunctionality < 0.5,] #0
sig_pathways_F_M_DEG_F <- sig_pathways_F_M_DEG_F %>% distinct(GeneMembers, .keep_all = T)
write.csv(sig_pathways_F_M_DEG_F, '2025_F_chrX_enriched_lipid_pathways.csv')

F_M_DEG_M_pr_out <- precRecall(annotation = female_male_DEG_M_reference, 
                    scores = female_male_DEG_M_pval,
                    scoreColumn = 2,
                    logTrans = T,
                    bigIsBetter = F,
                    iterations = 10000)

F_M_DEG_M_pathway_analysis <- as.data.frame(F_M_DEG_M_pr_out$results)

sig_pathways_F_M_DEG_M <- F_M_DEG_M_pathway_analysis[F_M_DEG_M_pathway_analysis$CorrectedPvalue < 0.05 & F_M_DEG_M_pathway_analysis$Multifunctionality < 0.5,] #0
sig_pathways_F_M_DEG_M <- sig_pathways_F_M_DEG_M %>% distinct(GeneMembers, .keep_all = T)
write.csv(sig_pathways_F_M_DEG_M, '2025_M_chrX_enriched_lipid_pathways.csv')



















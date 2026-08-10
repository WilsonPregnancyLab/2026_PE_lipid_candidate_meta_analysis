
#Differential Transcript Analysis using Swish

BiocManager::install("fishpond")
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("tximeta")
BiocManager::install("tximport")
BiocManager::install("org.Hs.eg.db")
library(fishpond) # version 1.6.0
library(org.Hs.eg.db) # version 3.23.1
library(tximport) # version 1.40.0
library(readr) # version 2.2.0
library(tximeta) # version 1.30.0
library(sva) # version 1.6.0
library(ggplot2) # version 3.5.1
library(gridExtra) # version 2.3
library(ggrepel) # version 0.9.5
library(SummarizedExperiment) # version 1.42.0
library(data.table) # version 1.18.4
library(tidyverse) # version 2.0.0


setwd("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/DESeq_rerun")

gtf19_file <- "./gencode.v37lift37.annotation.gtf"
gtf19_data <- import(gtf19_file)
gtf19_df <- as.data.frame(gtf19_data)

gtf_file <- "../genome_mapping/gencode.v48.primary_assembly.annotation.gtf"
gtf_data <- import(gtf_file)
gtf_df <- as.data.frame(gtf_data)
gtf_df_dups_rem <- gtf_df[!duplicated(gtf_df$transcript_id),]

lipid_gene_list <- read.csv("./lipid_genes_unique_P.csv")

setDT(lipid_gene_list)
setDT(gtf_df_dups_rem)
lipid_ensembl_gencode <- merge(x = lipid_gene_list, y = gtf_df_dups_rem[, c("gene_id", "gene_name", "seqnames", "gene_type", "transcript_id", "transcript_type", "transcript_name", "exon_number", "exon_id")], by.x = "gene_symbol", by.y = "gene_name", all.x = TRUE)
#lipid_ensembl <- named_read_counts[,c("Geneid", "gene_name")]
write.csv(lipid_ensembl_gencode, "lipid_ensembl_gencode.csv")

setwd ("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2026_Alternative_Splicing")
coldata <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/metadata/RNA_Lipid_Candidate_Metadata.csv")
coldata$exclude[coldata$predicted_fetal_sex == "F" & coldata$GSE_number == "GSE143953"] <- "exclude" #remove since not enough female only samples in this 
coldata$exclude[coldata$predicted_fetal_sex == "M" & coldata$GSE_number == "GSE148241"] <- "exclude"
coldata <- coldata[coldata$exclude != "exclude",]
coldata$files <- file.path("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2026_Alternative_Splicing", "quants", paste0("trimmed_", coldata$Run, "_quant"), "quant.sf")
all(file.exists(coldata$files))

##Read in quants w/ tximeta
library(SummarizedExperiment)

coldata$disease_group <- factor(coldata$disease_group, levels=c("Control", "PE"))
coldata$GSE_number <- factor(coldata$GSE_number)
#coldata$Run <- coldata$names
colnames(coldata)[colnames(coldata) == "Run"] <- "names"

se <- tximeta(coldata)
assayNames(se)
head(rownames(se))
y_se <- se
y_lip_se <- y_se[rownames(y_se) %in% lipid_ensembl_gencode$transcript_id,]

#Differential transcript expression - combined_fetal_sex
y <- scaleInfReps(y_se)
y <- labelKeep(y)
y <- y[mcols(y)$keep,]
set.seed(1)
y$combined_cov <- factor(paste0(y$GSE_number, "_", y$predicted_fetal_sex))
# colData(y)$GSE_number <- factor(colData(y)$GSE_number)


y <- swish(y, x="disease_group", cov="combined_cov")
y <- addIds(y, "SYMBOL", gene=TRUE)

res <- mcols(y)
head(res[order(res$qvalue), ])
table(res$qvalue < 0.05)
res_df <- data.frame(
    Transcript_ID = rownames(res),
    Tx_ID = res$tx_id,
    Gene_ID = res$gene_id,
    Log10_Mean = res$log10mean,
    Stat = res$stat,
    Log2_FoldChg = res$log2FC,
    P_Value = res$pvalue,
    Q_Value = res$qvalue,
    Locfdr = res$locfdr,
    Gene_Symbol = as.character(res$SYMBOL)
)
write.csv(res_df, "Swish_isoform_combined.csv")


# DTE_Swish <- function(df, covariate, file_suffix){
# y <- scaleInfReps(y_se)
# y <- labelKeep(y)
# y <- y[mcols(y)$keep,]
# set.seed(1)
# y <- swish(df, x="disease_group", cov= covariate)
# y <- addIds(y, "SYMBOL", gene=TRUE)

# res <- mcols(y)
# head(res[order(res$qvalue), ])
# table(res$qvalue < 0.05)
# res_df <- data.frame(
#     Transcript_ID = rownames(res),
#     Tx_ID = res$tx_id,
#     Gene_ID = res$gene_id,
#     Log10_Mean = res$log10mean,
#     Stat = res$stat,
#     Log2_FoldChg = res$log2FC,
#     P_Value = res$pvalue,
#     Q_Value = res$qvalue,
#     Locfdr = res$locfdr,
#     Gene_Symbol = as.character(res$SYMBOL)
# )
# filename <- paste0("Swish_isoform_", file_suffix, ".csv")
# write.csv(res_df, filename)
# }


# DTE_Swish(y,"combined_cov","combined")


# #Differential gene expression - combined_fetal_sex
# gse <- summarizeToGene(se)
# gy <- gse



# DGE_Swish <- function(df, covariate, file_suffix){

# gy <- scaleInfReps(gy)
# gy <- labelKeep(gy)
# gy <- gy[mcols(gy)$keep,]
# set.seed(1)
# gy <- swish(gy, x="disease_group", cov = "GSE_number")
# gy <- addIds(gy, "SYMBOL", gene=TRUE)
# table(mcols(gy)$qvalue < .05)

# res_gy <- mcols(gy)
# head(res_gy[order(res_gy$qvalue), ])
# table(res_gy$qvalue < 0.05)
# res_gy_df <- data.frame(
#     Gene_ID = res_gy$gene_id,
#     Stat = res_gy$stat,
#     Log2_FoldChg = res_gy$log2FC,
#     P_Value = res_gy$pvalue,
#     Q_Value = res_gy$qvalue,
#     Locfdr = res_gy$locfdr,
#     Gene_Symbol = as.character(res_gy$SYMBOL))

# }




#Differential transcript expression - combined_fetal_sex_lipid
y_lip <- scaleInfReps(y_lip_se)
y_lip <- labelKeep(y_lip)
y_lip <- y_lip[mcols(y_lip)$keep,]
set.seed(1)
y_lip$combined_cov <- factor(paste0(y_lip$GSE_number, "_", y_lip$predicted_fetal_sex))
# colData(y_lip)$GSE_number <- factor(colData(y_lip)$GSE_number)
y_lip <- swish(y_lip, x="disease_group", cov="combined_cov")
y_lip <- addIds(y_lip, "SYMBOL", gene=TRUE)

res <- mcols(y_lip)
head(res[order(res$qvalue), ])
table(res$qvalue < 0.05)
res_df <- data.frame(
    Transcript_ID = rownames(res),
    Tx_ID = res$tx_id,
    Gene_ID = res$gene_id,
    Log10_Mean = res$log10mean,
    Stat = res$stat,
    Log2_FoldChg = res$log2FC,
    P_Value = res$pvalue,
    Q_Value = res$qvalue,
    Locfdr = res$locfdr,
    Gene_Symbol = as.character(res$SYMBOL)
)
write.csv(res_df, "Swish_Lip_isoform_combined.csv")


#Differential transcript expression - fetal_female
y_se <- se
y_f_se <- y_se[, y_se$predicted_fetal_sex == "F"]
y_f <- scaleInfReps(y_f_se)
y_f <- labelKeep(y_f)
y_f <- y_f[mcols(y_f)$keep,]
set.seed(1)
colData(y_f)$GSE_number <- factor(colData(y_f)$GSE_number)
y_f <- swish(y_f, x="disease_group", cov="GSE_number")
y_f <- addIds(y_f, "SYMBOL", gene=TRUE)

res_f <- mcols(y_f)
head(res_f[order(res_f$qvalue), ])
table(res_f$qvalue < 0.05)
res_df_f <- data.frame(
    Transcript_ID = rownames(res_f),
    Tx_ID = res_f$tx_id,
    Gene_ID = res_f$gene_id,
    Log10_Mean = res_f$log10mean,
    Stat = res_f$stat,
    Log2_FoldChg = res_f$log2FC,
    P_Value = res_f$pvalue,
    Q_Value = res_f$qvalue,
    Locfdr = res_f$locfdr,
    Gene_Symbol = as.character(res_f$SYMBOL)
)
write.csv(res_df_f, "Swish_isoform_female.csv")

#Differential transcript expression - female_lipid
y_se <- se
y_lip_se <- y_se[rownames(y_se) %in% lipid_ensembl_gencode$transcript_id,]
y_f_lip_se <- y_lip_se[, y_lip_se$predicted_fetal_sex == "F"]
y_f_lip <- scaleInfReps(y_f_lip)
y_f_lip <- labelKeep(y_f_lip)
y_f_lip <- y_f_lip[mcols(y_f_lip)$keep,]
set.seed(1)
y_f_lip$combined_cov <- factor(paste0(y_f_lip$GSE_number, "_", y_f_lip$predicted_fetal_sex))
# colData(y_f_lip)$GSE_number <- factor(colData(y_f_lip)$GSE_number)
y_f_lip <- swish(y_f_lip, x="disease_group", cov="combined_cov")
y_f_lip <- addIds(y_f_lip, "SYMBOL", gene=TRUE)

res <- mcols(y_f_lip)
head(res[order(res$qvalue), ])
table(res$qvalue < 0.05)
res_df <- data.frame(
    Transcript_ID = rownames(res),
    Tx_ID = res$tx_id,
    Gene_ID = res$gene_id,
    Log10_Mean = res$log10mean,
    Stat = res$stat,
    Log2_FoldChg = res$log2FC,
    P_Value = res$pvalue,
    Q_Value = res$qvalue,
    Locfdr = res$locfdr,
    Gene_Symbol = as.character(res$SYMBOL)
)
write.csv(res_df, "Swish_Lip_isoform_female.csv")

#Differential transcript expression - fetal_male
y_se <- se
y_m_se <- y_se[, y_se$predicted_fetal_sex == "M"]
y_m <- scaleInfReps(y_m_se)
y_m <- labelKeep(y_m)
y_m <- y_m[mcols(y_m)$keep,]
set.seed(1)

colData(y_m)$GSE_number <- factor(colData(y_m)$GSE_number)
y_m <- swish(y_m, x="disease_group", cov="GSE_number")
y_m <- addIds(y_m, "SYMBOL", gene=TRUE)

res_m <- mcols(y_m)
head(res_m[order(res_m$qvalue), ])
table(res_m$qvalue < 0.05)
res_df_m <- data.frame(
    Transcript_ID = rownames(res_m),
    Tx_ID = res_m$tx_id,
    Gene_ID = res_m$gene_id,
    Log10_Mean = res_m$log10mean,
    Stat = res_m$stat,
    Log2_FoldChg = res_m$log2FC,
    P_Value = res_m$pvalue,
    Q_Value = res_m$qvalue,
    Locfdr = res_m$locfdr,
    Gene_Symbol = as.character(res_m$SYMBOL)
)
write.csv(res_df_m, "Swish_isoform_male.csv")


#Differential transcript expression - male_lipid
y_se <- se
y_lip_se <- y_se[rownames(y_se) %in% lipid_ensembl_gencode$transcript_id,]
y_m_lip_se <- y_lip_se[, y_lip_se$predicted_fetal_sex == "M"]
y_m_lip <- scaleInfReps(y_m_lip_se)
y_m_lip <- labelKeep(y_m_lip)
y_m_lip <- y_m_lip[mcols(y_m_lip)$keep,]
set.seed(1)
y_m_lip$combined_cov <- factor(paste0(y_m_lip$GSE_number, "_", y_m_lip$predicted_fetal_sex))
# colData(y_m_lip)$GSE_number <- factor(colData(y_m_lip)$GSE_number)
y_m_lip <- swish(y_m_lip, x="disease_group", cov="combined_cov")
y_m_lip <- addIds(y_m_lip, "SYMBOL", gene=TRUE)

res <- mcols(y_m_lip)
head(res[order(res$qvalue), ])
table(res$qvalue < 0.05)
res_df <- data.frame(
    Transcript_ID = rownames(res),
    Tx_ID = res$tx_id,
    Gene_ID = res$gene_id,
    Log10_Mean = res$log10mean,
    Stat = res$stat,
    Log2_FoldChg = res$log2FC,
    P_Value = res$pvalue,
    Q_Value = res$qvalue,
    Locfdr = res$locfdr,
    Gene_Symbol = as.character(res$SYMBOL)
)
write.csv(res_df, "Swish_Lip_isoform_male.csv")




#Use sex-specific rMATS
#Read alternative splicing results and combine to one table - combined_sex, F, M


rMATS_combine <- function(output_name, suffix){

dir_name <- paste0("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2026_Alternative_Splicing/alt_splic_analysis/", output_name)
setwd(dir_name)
SC_combined_JC <- read.table("./SE.MATS.JC.txt", header = TRUE, sep = "\t") #Skipped Exon (positive is more skipping in disease)
SC_combined_subset <- SC_combined_JC[,c("GeneID", "geneSymbol", "chr", "exonStart_0base", "exonEnd", "FDR", "IncLevelDifference")]
SC_combined_subset <- SC_combined_subset[SC_combined_subset$FDR <= 0.05 & (SC_combined_subset$IncLevelDifference <= -0.1 | SC_combined_subset$IncLevelDifference >= 0.1),]
colnames(SC_combined_subset)[colnames(SC_combined_subset) == "FDR"] <- "SC_FDR"
colnames(SC_combined_subset)[colnames(SC_combined_subset) == "IncLevelDifference"] <- "SC_IncLevelDifference"
colnames(SC_combined_subset)[colnames(SC_combined_subset) == "exonStart_0base"] <- "SC_exonStart_0base"
colnames(SC_combined_subset)[colnames(SC_combined_subset) == "exonEnd"] <- "SC_exonEnd"

RI_combined_JC <- read.table("./RI.MATS.JC.txt", header = TRUE, sep = "\t") #Retained Intron (positive is intron is more present in disease)
RI_combined_subset <- RI_combined_JC[,c("GeneID", "riExonStart_0base", "riExonEnd", "FDR", "IncLevelDifference")]
RI_combined_subset <- RI_combined_subset[RI_combined_subset$FDR <= 0.05 & (RI_combined_subset$IncLevelDifference <= -0.1 | RI_combined_subset$IncLevelDifference >= 0.1),]
colnames(RI_combined_subset)[colnames(RI_combined_subset) == "FDR"] <- "RI_FDR"
colnames(RI_combined_subset)[colnames(RI_combined_subset) == "IncLevelDifference"] <- "RI_IncLevelDifference"

MXE_combined_JC <- read.table("./MXE.MATS.JC.txt", header = TRUE, sep = "\t") #Mutually Exclusive Exons (positive is exon 1 is more often in controls, exon 2 is more often in PE)
MXE_combined_subset <- MXE_combined_JC[,c("GeneID", "X1stExonStart_0base", "X1stExonEnd", "X2ndExonStart_0base", "X2ndExonEnd", "FDR", "IncLevelDifference")]
MXE_combined_subset <- MXE_combined_subset[MXE_combined_subset$FDR <= 0.05 & (MXE_combined_subset$IncLevelDifference <= -0.1 | MXE_combined_subset$IncLevelDifference >= 0.1),]
colnames(MXE_combined_subset)[colnames(MXE_combined_subset) == "FDR"] <- "MXE_FDR"
colnames(MXE_combined_subset)[colnames(MXE_combined_subset) == "IncLevelDifference"] <- "MXE_IncLevelDifference"
colnames(MXE_combined_subset)[colnames(MXE_combined_subset) == "X1stExonStart_0base"] <- "MXE_X1stExonStart_0base"
colnames(MXE_combined_subset)[colnames(MXE_combined_subset) == "X1stExonEnd"] <- "MXE_X1stExonEnd"
colnames(MXE_combined_subset)[colnames(MXE_combined_subset) == "X2ndExonStart_0base"] <- "MXE_X2ndExonStart_0base"
colnames(MXE_combined_subset)[colnames(MXE_combined_subset) == "X2ndExonEnd"] <- "MXE_X2ndExonEnd"

A5SS_combined_JC <- read.table("./A5SS.MATS.JC.txt", header = TRUE, sep = "\t") #Alternative 5' Splice Site (positive is 5' splice site is shaved off on PE)
A5SS_combined_subset <- A5SS_combined_JC[,c("GeneID", "longExonStart_0base", "longExonEnd", "FDR", "IncLevelDifference")]
A5SS_combined_subset <- A5SS_combined_subset[A5SS_combined_subset$FDR <= 0.05 & (A5SS_combined_subset$IncLevelDifference <= -0.1 | A5SS_combined_subset$IncLevelDifference >= 0.1),]
colnames(A5SS_combined_subset)[colnames(A5SS_combined_subset) == "FDR"] <- "A5SS_FDR"
colnames(A5SS_combined_subset)[colnames(A5SS_combined_subset) == "IncLevelDifference"] <- "A5SS_IncLevelDifference"
colnames(A5SS_combined_subset)[colnames(A5SS_combined_subset) == "longExonStart_0base"] <- "A5SS_longExonStart_0base"
colnames(A5SS_combined_subset)[colnames(A5SS_combined_subset) == "longExonEnd"] <- "A5SS_longExonEnd"

A3SS_combined_JC <- read.table("./A3SS.MATS.JC.txt", header = TRUE, sep = "\t") #Alternative 3' Splicing Site (positive is 3' splice site is shaved off on PE)
A3SS_combined_subset <- A3SS_combined_JC[,c("GeneID", "longExonStart_0base", "longExonEnd", "FDR", "IncLevelDifference")]
A3SS_combined_subset <- A3SS_combined_subset[A3SS_combined_subset$FDR <= 0.05 & (A3SS_combined_subset$IncLevelDifference <= -0.1 | A3SS_combined_subset$IncLevelDifference >= 0.1),]
colnames(A3SS_combined_subset)[colnames(A3SS_combined_subset) == "FDR"] <- "A3SS_FDR"
colnames(A3SS_combined_subset)[colnames(A3SS_combined_subset) == "IncLevelDifference"] <- "A3SS_IncLevelDifference"
colnames(A3SS_combined_subset)[colnames(A3SS_combined_subset) == "longExonStart_0base"] <- "A3SS_longExonStart_0base"
colnames(A3SS_combined_subset)[colnames(A3SS_combined_subset) == "longExonEnd"] <- "A3SS_longExonEnd"

rMATS_list <- list(SC_combined_subset, RI_combined_subset, MXE_combined_subset, A5SS_combined_subset, A3SS_combined_subset)
rMATS_results <- rMATS_list %>% reduce(full_join, by = "GeneID")
rMATS_filename <- paste0("rMATS_results_", suffix, ".csv")

write.csv(rMATS_results, rMATS_filename)

}

rMATS_combine("output", "combined_sex")
rMATS_combine("output_F", "F")
rMATS_combine("output_M", "M")


# Part 1: NCBI Reference Subsetting
# the goal here is to add the NCBI annotations to the lipid genes (First, I'm going to do it for all of the lipid genes I tested - in whole_pop, M and F, then only for the significant ones in each group)

 if(!file.exists('Generic_human_ncbiIds_noParents.an.txt.gz')){
    system('wget https://gemma.msl.ubc.ca/annots/Generic_human_ncbiIds_noParents.an.txt.gz')}
NCBI <- read.table('Generic_human_ncbiIds_noParents.an.txt.gz', sep = '\t', header = T, quote = "")

#BiocManager::install("KEGGREST")
#BiocManager::install("clusterProfiler")
#devtools::install_github('PavlidisLab/ermineR')

library(KEGGREST)
library(clusterProfiler)
library(tidyverse)
library('org.Hs.eg.db') # get go terms from bioconductor 
library(ermineR) #version 1.0.3.9000
library(dplyr) #version 1.1.4
library(ggplot2) #version 3.5.1


# Whole Population Autosomes
WG_DESeq_auto_combined_sex <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_WG_rerun/RNA_WG_DESeq_results_autosomes_combined_sex.csv")
WG_DESeq_auto_combined_sex_pval <- WG_DESeq_auto_combined_sex[, c("gene_name", "pvalue")] #29251 genes
WG_DESeq_auto_combined_sex_pval <- WG_DESeq_auto_combined_sex_pval[!is.na(WG_DESeq_auto_combined_sex_pval$pvalue),] #29251
WG_DESeq_auto_combined_sex_pval <- WG_DESeq_auto_combined_sex_pval[!duplicated(WG_DESeq_auto_combined_sex_pval$gene_name),] #29141 genes investigated

WG_DESeq_auto_combined_sex_reference <- NCBI[NCBI$GeneSymbols %in% WG_DESeq_auto_combined_sex_pval$gene_name, ] #18375 genes

write.table(WG_DESeq_auto_combined_sex_pval, "WG_DESeq_auto_combined_sex_pval.txt", row.names = F, sep = '\t')
write.table(WG_DESeq_auto_combined_sex_reference, "WG_DESeq_auto_combined_sex_reference.txt", row.names = F, sep = '\t')

# Female Autosomes
RNA_seq_sex_stratified <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_WG_rerun/RNA_WG_DESeq_results_autosomes_IT_all.csv")
RNA_seq_autosomes_F <- RNA_seq_sex_stratified[RNA_seq_sex_stratified$comparison == "F_PEvsF_Cont", ]
# WG_DESeq_auto_F <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_WG_rerun/RNA_WG_DESeq_results_autosomes_F.csv")
WG_DESeq_auto_F_pval <- RNA_seq_autosomes_F[, c("gene_name", "pvalue")] 
WG_DESeq_auto_F_pval <- WG_DESeq_auto_F_pval[!is.na(WG_DESeq_auto_F_pval$pvalue),] 
WG_DESeq_auto_F_pval <- WG_DESeq_auto_F_pval[!duplicated(WG_DESeq_auto_F_pval$gene_name),] 

WG_DESeq_auto_F_reference <- NCBI[NCBI$GeneSymbols %in% WG_DESeq_auto_F_pval$gene_name, ] 

# Male Autosomes
RNA_seq_autosomes_M <- RNA_seq_sex_stratified[RNA_seq_sex_stratified$comparison == "M_PEvsM_Cont", ]
WG_DESeq_auto_M_pval <- RNA_seq_autosomes_M[, c("gene_name", "pvalue")] 
WG_DESeq_auto_M_pval <- WG_DESeq_auto_M_pval[!is.na(WG_DESeq_auto_M_pval$pvalue),] 
WG_DESeq_auto_M_pval <- WG_DESeq_auto_M_pval[!duplicated(WG_DESeq_auto_M_pval$gene_name),] 

WG_DESeq_auto_M_reference <- NCBI[NCBI$GeneSymbols %in% WG_DESeq_auto_M_pval$gene_name, ] 

# Female X Chromosomes
WG_DESeq_chrX_F <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_WG_rerun/RNA_WG_DESeq_results_chrX_F.csv")
WG_DESeq_chrX_F_pval <- WG_DESeq_chrX_F[, c("gene_name", "pvalue")] 
WG_DESeq_chrX_F_pval <- WG_DESeq_chrX_F_pval[!is.na(WG_DESeq_chrX_F_pval$pvalue),] 
WG_DESeq_chrX_F_pval <- WG_DESeq_chrX_F_pval[!duplicated(WG_DESeq_chrX_F_pval$gene_name),]

WG_DESeq_chrX_F_reference <- NCBI[NCBI$GeneSymbols %in% WG_DESeq_chrX_F_pval$gene_name, ] 

# Male X Chromosome
WG_DESeq_chrX_M <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_WG_rerun/RNA_WG_DESeq_results_chrX_M.csv")
WG_DESeq_chrX_M_pval <- WG_DESeq_chrX_M[, c("gene_name", "pvalue")] 
WG_DESeq_chrX_M_pval <- WG_DESeq_chrX_M_pval[!is.na(WG_DESeq_chrX_M_pval$pvalue),] 
WG_DESeq_chrX_M_pval <- WG_DESeq_chrX_M_pval[!duplicated(WG_DESeq_chrX_M_pval$gene_name),] 

WG_DESeq_chrX_M_reference <- NCBI[NCBI$GeneSymbols %in% WG_DESeq_chrX_M_pval$gene_name, ] 

# Part 2: Pathway Enrichment Analysis

Sys.setenv('JAVA_HOME' = '/usr/lib/jvm/java-21-openjdk/')

all_combinedpop_auto_pr_out <- precRecall(annotation = WG_DESeq_auto_combined_sex_reference, 
                    scores = WG_DESeq_auto_combined_sex_pval,
                    scoreColumn = 2,
                    logTrans = T,
                    bigIsBetter = F,
                    iterations = 10000)

all_combinedpop_auto_pathway_analysis <- as.data.frame(all_combinedpop_auto_pr_out$results)
sig_pathways_combinedpop_auto <- all_combinedpop_auto_pathway_analysis[all_combinedpop_auto_pathway_analysis$CorrectedPvalue < 0.05 & all_combinedpop_auto_pathway_analysis$Multifunctionality < 0.5,] #0
sig_pathways_combinedpop_auto <- sig_pathways_combinedpop_auto %>% distinct(GeneMembers, .keep_all = T)
write.csv(all_combinedpop_auto_pathway_analysis, '2025_combinedpop_auto_enriched_pathways.csv')

all_F_auto_pr_out <- precRecall(annotation = WG_DESeq_auto_F_reference, 
                    scores = WG_DESeq_auto_F_pval,
                    scoreColumn = 2,
                    logTrans = T,
                    bigIsBetter = F,
                    iterations = 10000)

all_F_auto_pathway_analysis <- as.data.frame(all_F_auto_pr_out$results)

sig_pathways_F_auto <- all_F_auto_pathway_analysis[all_F_auto_pathway_analysis$CorrectedPvalue < 0.05 & all_F_auto_pathway_analysis$Multifunctionality < 0.5,] #0
sig_pathways_F_auto <- sig_pathways_F_auto %>% distinct(GeneMembers, .keep_all = T)
write.csv(sig_pathways_F_auto, '2025_F_auto_enriched_pathways.csv')

all_M_auto_pr_out <- precRecall(annotation = WG_DESeq_auto_M_reference, 
                    scores = WG_DESeq_auto_M_pval,
                    scoreColumn = 2,
                    logTrans = T,
                    bigIsBetter = F,
                    iterations = 10000)

all_M_auto_pathway_analysis <- as.data.frame(all_M_auto_pr_out$results)

sig_pathways_M_auto <- all_M_auto_pathway_analysis[all_M_auto_pathway_analysis$CorrectedPvalue < 0.05 & all_M_auto_pathway_analysis$Multifunctionality < 0.5,] #0
sig_pathways_M_auto <- sig_pathways_M_auto %>% distinct(GeneMembers, .keep_all = T)
write.csv(sig_pathways_M_auto, '2025_M_auto_enriched_pathways.csv')

all_F_chrX_pr_out <- precRecall(annotation = WG_DESeq_chrX_F_reference, 
                    scores = WG_DESeq_chrX_F_pval,
                    scoreColumn = 2, 
                    logTrans = T,
                    bigIsBetter = F,
                    iterations = 10000)

all_F_chrX_pathway_analysis <- as.data.frame(all_F_chrX_pr_out$results)

sig_pathways_F_chrX <- all_F_chrX_pathway_analysis[all_F_chrX_pathway_analysis$CorrectedPvalue < 0.05 & all_F_chrX_pathway_analysis$Multifunctionality < 0.5,] #0
sig_pathways_F_chrX <- sig_pathways_F_chrX %>% distinct(GeneMembers, .keep_all = T)
write.csv(sig_pathways_F_chrX, '2025_F_chrX_enriched_pathways.csv')

all_M_chrX_pr_out <- precRecall(annotation = WG_DESeq_chrX_M_reference, 
                    scores = WG_DESeq_chrX_M_pval,
                    scoreColumn = 2,
                    logTrans = T,
                    bigIsBetter = F,
                    iterations = 10000)

all_M_chrX_pathway_analysis <- as.data.frame(all_M_chrX_pr_out$results)

sig_pathways_M_chrX <- all_M_chrX_pathway_analysis[all_M_chrX_pathway_analysis$CorrectedPvalue < 0.05 & all_M_chrX_pathway_analysis$Multifunctionality < 0.5,] #0
sig_pathways_M_chrX <- sig_pathways_M_chrX %>% distinct(GeneMembers, .keep_all = T)
write.csv(sig_pathways_M_chrX, '2025_M_chrX_enriched_pathways.csv')


#KEGG Pathway Enrichment


# Make KEGG Annotations
keggAnnots <- as.list(org.Hs.egPATH)
#keggAnnots <- keggAnnots %>% lapply(names)
keggAnnots %>% head

# Get Gene Symbols and Names (same as your GO code)
geneSymbols <- as.list(org.Hs.egSYMBOL) %>% unlist
geneName <- as.list(org.Hs.egGENENAME) %>% unlist

# Create your annotation object
kegg_annotation <- makeAnnotation(keggAnnots,
                            symbol = geneSymbols,
                            name = geneName,
                            output = NULL, 
                            return = TRUE)

names(kegg_annotation)[names(kegg_annotation) == "GOTerms"] <- "KEGG_Pathway"

# Whole Population Autosomes

KEGG_WG_DESeq_auto_combined_sex_reference <- kegg_annotation[kegg_annotation$GeneSymbols %in% WG_DESeq_auto_combined_sex_pval$gene_name, ] #18375 genes

# Female Autosomes

KEGG_WG_DESeq_auto_F_reference <- kegg_annotation[kegg_annotation$GeneSymbols %in% WG_DESeq_auto_F_pval$gene_name, ] 

# Male Autosomes

KEGG_WG_DESeq_auto_M_reference <- kegg_annotation[kegg_annotation$GeneSymbols %in% WG_DESeq_auto_M_pval$gene_name, ] 

# Female X Chromosomes

KEGG_WG_DESeq_chrX_F_reference <- kegg_annotation[kegg_annotation$GeneSymbols %in% WG_DESeq_chrX_F_pval$gene_name, ] 

# Male X Chromosome

KEGG_WG_DESeq_chrX_M_reference <- kegg_annotation[kegg_annotation$GeneSymbols %in% WG_DESeq_chrX_M_pval$gene_name, ] 

# Part 2: Pathway Enrichment Analysis

Sys.setenv('JAVA_HOME' = '/usr/lib/jvm/java-21-openjdk/')

KEGG_all_combinedpop_auto_pr_out <- precRecall(annotation = KEGG_WG_DESeq_auto_combined_sex_reference, 
                    scores = WG_DESeq_auto_combined_sex_pval,
                    scoreColumn = 2,
                    logTrans = T,
                    bigIsBetter = F,
                    iterations = 10000)

KEGG_all_combinedpop_auto_pathway_analysis <- as.data.frame(KEGG_all_combinedpop_auto_pr_out$results)
KEGG_sig_pathways_combinedpop_auto <- KEGG_all_combinedpop_auto_pathway_analysis[KEGG_all_combinedpop_auto_pathway_analysis$CorrectedPvalue < 0.05 & KEGG_all_combinedpop_auto_pathway_analysis$Multifunctionality < 0.5,] #0
KEGG_sig_pathways_combinedpop_auto <- KEGG_sig_pathways_combinedpop_auto %>% distinct(GeneMembers, .keep_all = T)
write.csv(KEGG_all_combinedpop_auto_pathway_analysis, 'KEGG_2025_combinedpop_auto_enriched_pathways.csv')

KEGG_all_F_auto_pr_out <- precRecall(annotation = KEGG_WG_DESeq_auto_F_reference, 
                    scores = WG_DESeq_auto_F_pval,
                    scoreColumn = 2,
                    logTrans = T,
                    bigIsBetter = F,
                    iterations = 10000)

KEGG_all_F_auto_pathway_analysis <- as.data.frame(KEGG_all_F_auto_pr_out$results)

KEGG_sig_pathways_F_auto <- KEGG_all_F_auto_pathway_analysis[KEGG_all_F_auto_pathway_analysis$CorrectedPvalue < 0.05 & KEGG_all_F_auto_pathway_analysis$Multifunctionality < 0.5,] #0
KEGG_sig_pathways_F_auto <- KEGG_sig_pathways_F_auto %>% distinct(GeneMembers, .keep_all = T)
write.csv(KEGG_sig_pathways_F_auto, 'KEGG_2025_F_auto_enriched_pathways.csv')

KEGG_all_M_auto_pr_out <- precRecall(annotation = KEGG_WG_DESeq_auto_M_reference, 
                    scores = WG_DESeq_auto_M_pval,
                    scoreColumn = 2,
                    logTrans = T,
                    bigIsBetter = F,
                    iterations = 10000)

KEGG_all_M_auto_pathway_analysis <- as.data.frame(KEGG_all_M_auto_pr_out$results)

KEGG_sig_pathways_M_auto <- KEGG_all_M_auto_pathway_analysis[KEGG_all_M_auto_pathway_analysis$CorrectedPvalue < 0.05 & KEGG_all_M_auto_pathway_analysis$Multifunctionality < 0.5,] #0
KEGG_sig_pathways_M_auto <- KEGG_sig_pathways_M_auto %>% distinct(GeneMembers, .keep_all = T)
write.csv(KEGG_sig_pathways_M_auto, 'KEGG_2025_M_auto_enriched_pathways_KEGG.csv')

KEGG_all_F_chrX_pr_out <- precRecall(annotation = KEGG_WG_DESeq_chrX_F_reference, 
                    scores = WG_DESeq_chrX_F_pval,
                    scoreColumn = 2, 
                    logTrans = T,
                    bigIsBetter = F,
                    iterations = 10000)

KEGG_all_F_chrX_pathway_analysis <- as.data.frame(KEGG_all_F_chrX_pr_out$results)

KEGG_sig_pathways_F_chrX <- KEGG_all_F_chrX_pathway_analysis[KEGG_all_F_chrX_pathway_analysis$CorrectedPvalue < 0.05 & KEGG_all_F_chrX_pathway_analysis$Multifunctionality < 0.5,] #0
KEGG_sig_pathways_F_chrX <- KEGG_sig_pathways_F_chrX %>% distinct(GeneMembers, .keep_all = T)
write.csv(KEGG_sig_pathways_F_chrX, 'KEGG_2025_F_chrX_enriched_pathways.csv')

KEGG_all_M_chrX_pr_out <- precRecall(annotation = KEGG_WG_DESeq_chrX_M_reference, 
                    scores = WG_DESeq_chrX_M_pval,
                    scoreColumn = 2,
                    logTrans = T,
                    bigIsBetter = F,
                    iterations = 10000)

KEGG_all_M_chrX_pathway_analysis <- as.data.frame(KEGG_all_M_chrX_pr_out$results)

KEGG_sig_pathways_M_chrX <- KEGG_all_M_chrX_pathway_analysis[KEGG_all_M_chrX_pathway_analysis$CorrectedPvalue < 0.05 & KEGG_all_M_chrX_pathway_analysis$Multifunctionality < 0.5,] #0
KEGG_sig_pathways_M_chrX <- KEGG_sig_pathways_M_chrX %>% distinct(GeneMembers, .keep_all = T)
write.csv(KEGG_sig_pathways_M_chrX, 'KEGG_2025_M_chrX_enriched_pathways.csv')


















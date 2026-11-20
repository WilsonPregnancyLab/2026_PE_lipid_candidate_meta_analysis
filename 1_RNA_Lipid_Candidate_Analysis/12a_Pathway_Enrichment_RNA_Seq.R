devtools::install_github('PavlidisLab/ermineR')
library(ermineR) #version 1.0.3.9000
library(dplyr) #version 1.1.4
library(ggplot2) #version 3.5.1
library(org.Hs.eg.db)

# Part 1: NCBI Reference Subsetting
# the goal here is to add the NCBI annotations to the lipid genes (First, I'm going to do it for all of the lipid genes I tested - in whole_pop, M and F, then only for the significant ones in each group)

 if(!file.exists('Generic_human_ncbiIds_noParents.an.txt.gz')){
    system('wget https://gemma.msl.ubc.ca/annots/Generic_human_ncbiIds_noParents.an.txt.gz')}
NCBI <- read.table('Generic_human_ncbiIds_noParents.an.txt.gz', sep = '\t', header = T, quote = "")

# Whole Population Autosomes
Lip_DESeq_auto_combined_sex <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/DESeq_rerun/RNA_Lip_DESeq_results_autosomes_combined_sex.csv")
Lip_DESeq_auto_combined_sex_pval <- Lip_DESeq_auto_combined_sex[, c("gene_symbol", "pvalue")] #4750 genes
Lip_DESeq_auto_combined_sex_pval <- Lip_DESeq_auto_combined_sex_pval[!is.na(Lip_DESeq_auto_combined_sex_pval$pvalue),] #4944
Lip_DESeq_auto_combined_sex_pval <- Lip_DESeq_auto_combined_sex_pval[!duplicated(Lip_DESeq_auto_combined_sex_pval$gene_symbol),] #4750 lipid genes investigated

Lip_DESeq_auto_combined_sex_reference <- NCBI[NCBI$GeneSymbols %in% Lip_DESeq_auto_combined_sex_pval$gene_symbol, ] #4739 genes

write.table(Lip_DESeq_auto_combined_sex_pval, "Lip_DESeq_auto_combined_sex_pval.txt", row.names = F, sep = '\t')
write.table(Lip_DESeq_auto_combined_sex_reference, "Lip_DESeq_auto_combined_sex_reference.txt", row.names = F, sep = '\t')


# Female Autosomes
RNA_seq_sex_stratified <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_WG_rerun/RNA_WG_DESeq_results_autosomes_IT_all.csv")
RNA_seq_autosomes_F <- RNA_seq_sex_stratified[RNA_seq_sex_stratified$comparison == "F_PEvsF_Cont", ]
Lip_DESeq_auto_F_pval <- RNA_seq_autosomes_F[, c("gene_name", "pvalue")] # 4710 genes
Lip_DESeq_auto_F_pval <- Lip_DESeq_auto_F_pval[!is.na(Lip_DESeq_auto_F_pval$pvalue),] #4710
Lip_DESeq_auto_F_pval <- Lip_DESeq_auto_F_pval[!duplicated(Lip_DESeq_auto_F_pval$gene_name),] # 4708 lipid genes investigated

Lip_DESeq_auto_F_reference <- NCBI[NCBI$GeneSymbols %in% Lip_DESeq_auto_F_pval$gene_name, ] #4698 reference genes

# Male Autosomes
RNA_seq_autosomes_M <- RNA_seq_sex_stratified[RNA_seq_sex_stratified$comparison == "M_PEvsM_Cont", ]
Lip_DESeq_auto_M_pval <- RNA_seq_autosomes_M[, c("gene_name", "pvalue")] # 4614 probes
Lip_DESeq_auto_M_pval <- Lip_DESeq_auto_M_pval[!is.na(Lip_DESeq_auto_M_pval$pvalue),] #4858
Lip_DESeq_auto_M_pval <- Lip_DESeq_auto_M_pval[!duplicated(Lip_DESeq_auto_M_pval$gene_name),] # 4613 lipid genes investigated

Lip_DESeq_auto_M_reference <- NCBI[NCBI$GeneSymbols %in% Lip_DESeq_auto_M_pval$gene_name, ] #4607 genes

# Female X Chromosomes
Lip_DESeq_chrX_F <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/DESeq_rerun/RNA_Lip_DESeq_results_chrX_F.csv")
Lip_DESeq_chrX_F_pval <- Lip_DESeq_chrX_F[, c("gene_symbol", "pvalue")] # 149 genes
Lip_DESeq_chrX_F_pval <- Lip_DESeq_chrX_F_pval[!is.na(Lip_DESeq_chrX_F_pval$pvalue),] #149
Lip_DESeq_chrX_F_pval <- Lip_DESeq_chrX_F_pval[!duplicated(Lip_DESeq_chrX_F_pval$gene_symbol),] # 149 lipid genes investigated

Lip_DESeq_chrX_F_reference <- NCBI[NCBI$GeneSymbols %in% Lip_DESeq_chrX_F_pval$gene_symbol, ] # 149 reference genes

# Male X Chromosome
Lip_DESeq_chrX_M <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/DESeq_rerun/RNA_Lip_DESeq_results_chrX_M.csv")
Lip_DESeq_chrX_M_pval <- Lip_DESeq_chrX_M[, c("gene_symbol", "pvalue")] # 147 probes
Lip_DESeq_chrX_M_pval <- Lip_DESeq_chrX_M_pval[!is.na(Lip_DESeq_chrX_M_pval$pvalue),] #147
Lip_DESeq_chrX_M_pval <- Lip_DESeq_chrX_M_pval[!duplicated(Lip_DESeq_chrX_M_pval$gene_symbol),] # 147 lipid genes investigated

Lip_DESeq_chrX_M_reference <- NCBI[NCBI$GeneSymbols %in% Lip_DESeq_chrX_M_pval$gene_symbol, ] #147 genes

# Part 2: Pathway Enrichment Analysis

Sys.setenv('JAVA_HOME' = '/usr/lib/jvm/java-21-openjdk/')

Lip_DESeq_auto_combined_sex_pval <- read.table("Lip_DESeq_auto_combined_sex_pval.txt")
Lip_DESeq_auto_combined_sex_reference <- read.table("Lip_DESeq_auto_combined_sex_reference.txt")
Lip_DESeq_auto_F_pval <- read.csv("Lip_DESeq_auto_F_pval.csv")
Lip_DESeq_auto_M_pval <- read.csv("Lip_DESeq_auto_M_pval.csv")
Lip_DESeq_chrX_F_pval <- read.csv("Lip_DESeq_chrX_F_pval.csv")
Lip_DESeq_chrX_M_pval <- read.csv("Lip_DESeq_chrX_M_pval.csv")


all_combinedpop_auto_pr_out <- precRecall(annotation = Lip_DESeq_auto_combined_sex_reference, 
                    scores = Lip_DESeq_auto_combined_sex_pval,
                    scoreColumn = 2,
                    logTrans = T,
                    bigIsBetter = F,
                    iterations = 10000)

all_combinedpop_auto_pathway_analysis <- as.data.frame(all_combinedpop_auto_pr_out$results)
sig_pathways_combinedpop_auto <- all_combinedpop_auto_pathway_analysis[all_combinedpop_auto_pathway_analysis$CorrectedPvalue < 0.05 & all_combinedpop_auto_pathway_analysis$Multifunctionality < 0.5,] #0
sig_pathways_combinedpop_auto <- sig_pathways_combinedpop_auto %>% distinct(GeneMembers, .keep_all = T)
write.csv(all_combinedpop_auto_pathway_analysis, '2025_combinedpop_auto_enriched_lipid_pathways.csv')

all_F_auto_pr_out <- precRecall(annotation = Lip_DESeq_auto_F_reference, 
                    scores = Lip_DESeq_auto_F_pval,
                    scoreColumn = 2,
                    logTrans = T,
                    bigIsBetter = F,
                    iterations = 10000)

all_F_auto_pathway_analysis <- as.data.frame(all_F_auto_pr_out$results)

sig_pathways_F_auto <- all_F_auto_pathway_analysis[all_F_auto_pathway_analysis$CorrectedPvalue < 0.05 & all_F_auto_pathway_analysis$Multifunctionality < 0.5,] #0
sig_pathways_F_auto <- sig_pathways_F_auto %>% distinct(GeneMembers, .keep_all = T)
write.csv(sig_pathways_F_auto, '2025_F_auto_enriched_lipid_pathways.csv')

all_M_auto_pr_out <- precRecall(annotation = Lip_DESeq_auto_M_reference, 
                    scores = Lip_DESeq_auto_M_pval,
                    scoreColumn = 2,
                    logTrans = T,
                    bigIsBetter = F,
                    iterations = 10000)

all_M_auto_pathway_analysis <- as.data.frame(all_M_auto_pr_out$results)

sig_pathways_M_auto <- all_M_auto_pathway_analysis[all_M_auto_pathway_analysis$CorrectedPvalue < 0.05 & all_M_auto_pathway_analysis$Multifunctionality < 0.5,] #0
sig_pathways_M_auto <- sig_pathways_M_auto %>% distinct(GeneMembers, .keep_all = T)
write.csv(sig_pathways_M_auto, '2025_M_auto_enriched_lipid_pathways.csv')

all_F_chrX_pr_out <- precRecall(annotation = Lip_DESeq_chrX_F_reference, 
                    scores = Lip_DESeq_chrX_F_pval,
                    scoreColumn = 2, 
                    logTrans = T,
                    bigIsBetter = F,
                    iterations = 10000)

all_F_chrX_pathway_analysis <- as.data.frame(all_F_chrX_pr_out$results)

sig_pathways_F_chrX <- all_F_chrX_pathway_analysis[all_F_chrX_pathway_analysis$CorrectedPvalue < 0.05 & all_F_chrX_pathway_analysis$Multifunctionality < 0.5,] #0
sig_pathways_F_chrX <- sig_pathways_F_chrX %>% distinct(GeneMembers, .keep_all = T)
write.csv(sig_pathways_F_chrX, '2025_F_chrX_enriched_lipid_pathways.csv')

all_M_chrX_pr_out <- precRecall(annotation = Lip_DESeq_chrX_M_reference, 
                    scores = Lip_DESeq_chrX_M_pval,
                    scoreColumn = 2,
                    logTrans = T,
                    bigIsBetter = F,
                    iterations = 10000)

all_M_chrX_pathway_analysis <- as.data.frame(all_M_chrX_pr_out$results)

sig_pathways_M_chrX <- all_M_chrX_pathway_analysis[all_M_chrX_pathway_analysis$CorrectedPvalue < 0.05 & all_M_chrX_pathway_analysis$Multifunctionality < 0.5,] #0
sig_pathways_M_chrX <- sig_pathways_M_chrX %>% distinct(GeneMembers, .keep_all = T)
write.csv(sig_pathways_M_chrX, '2025_M_chrX_enriched_lipid_pathways.csv')





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

KEGG_Lip_DESeq_auto_combined_sex_reference <- kegg_annotation[kegg_annotation$GeneSymbols %in% Lip_DESeq_auto_combined_sex_pval$gene_symbol, ] #18375 genes

# Female Autosomes

KEGG_Lip_DESeq_auto_F_reference <- kegg_annotation[kegg_annotation$GeneSymbols %in% Lip_DESeq_auto_F_pval$gene_name, ] 

# Male Autosomes

KEGG_Lip_DESeq_auto_M_reference <- kegg_annotation[kegg_annotation$GeneSymbols %in% Lip_DESeq_auto_M_pval$gene_name, ] 

# Female X Chromosomes

KEGG_Lip_DESeq_chrX_F_reference <- kegg_annotation[kegg_annotation$GeneSymbols %in% Lip_DESeq_chrX_F_pval$gene_symbol, ] 

# Male X Chromosome

KEGG_Lip_DESeq_chrX_M_reference <- kegg_annotation[kegg_annotation$GeneSymbols %in% Lip_DESeq_chrX_M_pval$gene_symbol, ] 

# Part 2: Pathway Enrichment Analysis

Sys.setenv('JAVA_HOME' = '/usr/lib/jvm/java-21-openjdk/')

KEGG_all_combinedpop_auto_pr_out <- precRecall(annotation = KEGG_Lip_DESeq_auto_combined_sex_reference, 
                    scores = Lip_DESeq_auto_combined_sex_pval,
                    scoreColumn = 2,
                    logTrans = T,
                    bigIsBetter = F,
                    iterations = 10000)

KEGG_all_combinedpop_auto_pathway_analysis <- as.data.frame(KEGG_all_combinedpop_auto_pr_out$results)
KEGG_sig_pathways_combinedpop_auto <- KEGG_all_combinedpop_auto_pathway_analysis[KEGG_all_combinedpop_auto_pathway_analysis$CorrectedPvalue < 0.05 & KEGG_all_combinedpop_auto_pathway_analysis$Multifunctionality < 0.5,] #0
KEGG_sig_pathways_combinedpop_auto <- KEGG_sig_pathways_combinedpop_auto %>% distinct(GeneMembers, .keep_all = T)
write.csv(KEGG_all_combinedpop_auto_pathway_analysis, 'KEGG_2025_combinedpop_auto_enriched_pathways.csv')

KEGG_all_F_auto_pr_out <- precRecall(annotation = KEGG_Lip_DESeq_auto_F_reference, 
                    scores = Lip_DESeq_auto_F_pval,
                    scoreColumn = 2,
                    logTrans = T,
                    bigIsBetter = F,
                    iterations = 10000)

KEGG_all_F_auto_pathway_analysis <- as.data.frame(KEGG_all_F_auto_pr_out$results)
KEGG_sig_pathways_F_auto <- KEGG_all_F_auto_pathway_analysis[KEGG_all_F_auto_pathway_analysis$CorrectedPvalue < 0.05 & KEGG_all_F_auto_pathway_analysis$Multifunctionality < 0.5,] #0
KEGG_sig_pathways_F_auto <- KEGG_sig_pathways_F_auto %>% distinct(GeneMembers, .keep_all = T)
write.csv(KEGG_sig_pathways_F_auto, 'KEGG_2025_F_auto_enriched_pathways.csv')

KEGG_all_M_auto_pr_out <- precRecall(annotation = KEGG_Lip_DESeq_auto_M_reference, 
                    scores = Lip_DESeq_auto_M_pval,
                    scoreColumn = 2,
                    logTrans = T,
                    bigIsBetter = F,
                    iterations = 10000)

KEGG_all_M_auto_pathway_analysis <- as.data.frame(KEGG_all_M_auto_pr_out$results)

KEGG_sig_pathways_M_auto <- KEGG_all_M_auto_pathway_analysis[KEGG_all_M_auto_pathway_analysis$CorrectedPvalue < 0.05 & KEGG_all_M_auto_pathway_analysis$Multifunctionality < 0.5,] #0
KEGG_sig_pathways_M_auto <- KEGG_sig_pathways_M_auto %>% distinct(GeneMembers, .keep_all = T)
write.csv(KEGG_sig_pathways_M_auto, 'KEGG_2025_M_auto_enriched_pathways_KEGG.csv')

KEGG_all_F_chrX_pr_out <- precRecall(annotation = KEGG_Lip_DESeq_chrX_F_reference, 
                    scores = Lip_DESeq_chrX_F_pval,
                    scoreColumn = 2, 
                    logTrans = T,
                    bigIsBetter = F,
                    iterations = 10000)

KEGG_all_F_chrX_pathway_analysis <- as.data.frame(KEGG_all_F_chrX_pr_out$results)

KEGG_sig_pathways_F_chrX <- KEGG_all_F_chrX_pathway_analysis[KEGG_all_F_chrX_pathway_analysis$CorrectedPvalue < 0.05 & KEGG_all_F_chrX_pathway_analysis$Multifunctionality < 0.5,] #0
KEGG_sig_pathways_F_chrX <- KEGG_sig_pathways_F_chrX %>% distinct(GeneMembers, .keep_all = T)
write.csv(KEGG_sig_pathways_F_chrX, 'KEGG_2025_F_chrX_enriched_pathways.csv')

KEGG_all_M_chrX_pr_out <- precRecall(annotation = KEGG_Lip_DESeq_chrX_M_reference, 
                    scores = Lip_DESeq_chrX_M_pval,
                    scoreColumn = 2,
                    logTrans = T,
                    bigIsBetter = F,
                    iterations = 10000)

KEGG_all_M_chrX_pathway_analysis <- as.data.frame(KEGG_all_M_chrX_pr_out$results)

KEGG_sig_pathways_M_chrX <- KEGG_all_M_chrX_pathway_analysis[KEGG_all_M_chrX_pathway_analysis$CorrectedPvalue < 0.05 & KEGG_all_M_chrX_pathway_analysis$Multifunctionality < 0.5,] #0
KEGG_sig_pathways_M_chrX <- KEGG_sig_pathways_M_chrX %>% distinct(GeneMembers, .keep_all = T)
write.csv(KEGG_sig_pathways_M_chrX, 'KEGG_2025_M_chrX_enriched_pathways.csv')







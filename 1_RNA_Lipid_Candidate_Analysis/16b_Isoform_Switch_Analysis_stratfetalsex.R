BiocManager::install("IsoformSwitchAnalyzeR")
library(IsoformSwitchAnalyzeR) #version 2.12.0
library(rtracklayer) #version 1.72.0
library(tidyverse) #version 2.0.0
library(dplyr) #version 1.2.1
library(purrr) #version 1.2.2
library(patchwork) #version 1.3.2
library(ggplot2) #version 4.0.3
library(scales) #version 1.4.0, Needed for percentage formatting on the y-axis


# like the PTC analysis, alternative splicing analysis, prediction of consequences

#Importing Data
salmonQuant <- importIsoformExpression(
    parentDir = "/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2026_Alternative_Splicing/quants"
)

head(salmonQuant$abundance, 2)
head(salmonQuant$counts, 2) #395452 2

#subset the transcript IDs that belong ONLY to the Primary Assembly GTF
gtf_data <- import("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/genome_mapping/gencode.v48.primary_assembly.annotation.gtf")
primary_transcripts <- unique(gtf_data$transcript_id)
primary_transcripts <- na.omit(primary_transcripts) #387954
salmonQuant_f <- salmonQuant
salmonQuant_f$counts    <- salmonQuant$counts[(salmonQuant$counts)$isoform_id %in% primary_transcripts, ]
salmonQuant_f$abundance <- salmonQuant$abundance[(salmonQuant$abundance)$isoform_id %in% primary_transcripts, ]

metadata <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/metadata/RNA_Lipid_Candidate_Metadata.csv")
metadata <- metadata[metadata$exclude != "exclude",]
metadata$files <- file.path("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2026_Alternative_Splicing", "quants", paste0("trimmed_", metadata$Run, "_quant"), "quant.sf")
metadata$sampleID <-  paste0("trimmed_", metadata$Run, "_quant")
metadata$cond_sex <- paste0(metadata$disease_group, "_", metadata$predicted_fetal_sex)

columns <- c("sampleID", "GSE_number", "cond_sex")
myDesign <- metadata[, colnames(metadata) %in% columns]
myDesign_ordered <- myDesign[, c("sampleID", "cond_sex", "GSE_number")]
colnames(myDesign_ordered)[colnames(myDesign_ordered) == "cond_sex"] <- "condition"
myDesign_ordered$sampleID   <- as.character(myDesign_ordered$sampleID)
myDesign_ordered$condition  <- as.character(myDesign_ordered$condition)
myDesign_ordered$GSE_number <- as.character(myDesign_ordered$GSE_number)
# myDesign_ordered$predicted_fetal_sex <- as.character(myDesign_ordered$predicted_fetal_sex)

# all(file.exists(myDesign$files))

aSwitchList_fet_sex <- importRdata(
    isoformCountMatrix    = salmonQuant_f$counts,
    isoformRepExpression  = salmonQuant_f$abundance,
    designMatrix          = myDesign_ordered,
    isoformExonAnnoation = "/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/genome_mapping/gencode.v48.primary_assembly.annotation.gtf",
    isoformNtFasta        = "/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2026_Alternative_Splicing/gencode.v48.transcripts.fa",
    fixStringTieAnnotationProblem = TRUE,
    removeNonConvensionalChr = TRUE,
    showProgress = TRUE,
    estimateDifferentialGeneRange=FALSE)


summary(aSwitchList)

#Identify Isoform Switches - Whole Genome

dir.create("./Isoform_DEXSeq_fet_sex")
aSwitchList_fet_sex <- preFilter(aSwitchList_fet_sex,
                isoCount = 10,
                min.Count.prop = 0.7,
                IFcutoff = 0.1,
                min.IF.prop = 0.5,
                alpha = 0.05,
                dIFcutoff = 0.1) 

exampleSwitchListAnalyzed_fet_sex <- isoformSwitchTestDEXSeq(
    switchAnalyzeRlist   = aSwitchList_fet_sex,
    reduceToSwitchingGenes = FALSE
)

#Summary of isoform switches
Isoform_Switch_Summary_fet_sex <- extractSwitchSummary(exampleSwitchListAnalyzed_fet_sex)
write.csv(Isoform_Switch_Summary_fet_sex, "Isoform_Switch_Summary_fet_sex_WG.csv")

exampleSwitchListAnalyzed_fet_sex <- extractSequence(
    exampleSwitchListAnalyzed_fet_sex, 
    pathToOutput = "./Isoform_DEXSeq_fet_sex"
)

#In browser paste isoformSwitchAnalyzeR_isoform_nt.fasta file into https://cpc2.gao-lab.org/ to get information about coding potential (CPC2 - Coding Potential Calculator 2)

#In browser paste isoformSwitchAnalyzeR_isoform_AA.fasta file into https://www.ebi.ac.uk/jdispatcher/pfa/pfamscan to get information about protein families impacted

#In browser paste isoformSwitchAnalyzeR_isoform_AA.fasta file into https://services.healthtech.dtu.dk/services/SignalP-6.0/ to predict signal peptides and their cleavage sites

#In browser paste isoformSwitchAnalyzeR_isoform_AA.fasta file into https://services.healthtech.dtu.dk/services/SignalP-6.0/ to predict signal peptides and their cleavage sites


dir.create("./Isoform_P2_fet_sex")

##Copy all those results .txt files into "./Isoform_P2" directory

#Analysis of Alternative Splicing and functional consequences

exampleSwitchListAnalyzed_fet_sex <- isoformSwitchAnalysisPart2(
  switchAnalyzeRlist        = exampleSwitchListAnalyzed_fet_sex, 
  n                         = 50,    # if plotting was enabled, it would only output the top 10 switches
  removeNoncodinORFs        = FALSE,
  pathToCPC2resultFile      = "./Isoform_DEXSeq_fet_sex/cpc2_result.txt", #adds column identifying if transcript is protein- or non-coding
  pathToPFAMresultFile      = "./Isoform_DEXSeq_fet_sex/pfam_results.txt", #identifies functional and structural domains located on exons
  pathToIUPred2AresultFile  = "./Isoform_DEXSeq_fet_sex/iupred2a_result.txt", 
  pathToSignalPresultFile   = "./Isoform_DEXSeq_fet_sex/signalP_results.txt", #marks disordered regions (involved in cell signaling)
  outputPlots               = TRUE)

exampleSwitchListAnalyzed_fet_sex <- analyzeAlternativeSplicing(
  switchAnalyzeRlist = exampleSwitchListAnalyzed_fet_sex,
  quiet = FALSE
)

saveRDS(exampleSwitchListAnalyzed_fet_sex, file = "./Isoform_DEXSeq/exampleSwitchListAnalyzed_fet_sex.rds")

# write.csv(exampleSwitchListAnalyzed_fet_sex$isoformFeatures,"PE_isoform_features_fet_sex.csv", row.names = FALSE)

#Only females had sig_diff isoforms and Switches
 WG_female_SwitchList <- subsetSwitchAnalyzeRlist(
  switchAnalyzeRlist = exampleSwitchListAnalyzed_fet_sex,
  subset = exampleSwitchListAnalyzed_fet_sex$isoformFeatures$condition_1 == "Control_F" &
           exampleSwitchListAnalyzed_fet_sex$isoformFeatures$condition_2 == "PE_F"
)

 WG_male_SwitchList <- subsetSwitchAnalyzeRlist(
  switchAnalyzeRlist = exampleSwitchListAnalyzed_fet_sex,
  subset = exampleSwitchListAnalyzed_fet_sex$isoformFeatures$condition_1 == "Control_M" &
           exampleSwitchListAnalyzed_fet_sex$isoformFeatures$condition_2 == "PE_M"
)

dir.create("./WG_F_Isoform_Tables")
setwd("./WG_F_Isoform_Tables")
Isoform_Table(WG_female_SwitchList, "WG_F", RNA_Seq_isoform_female, RNA_WG_DESeq_results_autosomes_F, rMATS_results_WG)

setwd("../")
dir.create("./WG_M_Isoform_Tables")
setwd("./WG_M_Isoform_Tables")
Isoform_Table(WG_male_SwitchList, "WG_M", RNA_Seq_isoform_male, RNA_WG_DESeq_results_autosomes_M, rMATS_results_WG)


#Plots

png("./Isoform_Splicing_Summary_fet_sex.png", height = 9, width = 15, units = "in", res = 500)
extractSplicingSummary(
    exampleSwitchListAnalyzed_fet_sex,
    asFractionTotal = FALSE,
    plotGenes=FALSE)
dev.off()

png("./Isoform_Splicing_Enrichment_fet_sex.png", height = 9, width = 15, units = "in", res = 500)
extractSplicingEnrichment(
    exampleSwitchListAnalyzed_fet_sex,
    returnResult = TRUE # if TRUE returns a data.frame with the summary statistics
)
dev.off() # NO SPLICING ENRICHMENT FOR MY DATASET



png("./Isoform_Switching_Summary_fet_sex.png", height = 9, width = 15, units = "in", res = 500)
extractTopSwitches(exampleSwitchListAnalyzed_fet_sex, filterForConsequences = TRUE)
dev.off()


#Consequence of Isoform Switches
png("./Isoform_Consequence_Summary_fet_sex.png", height = 9, width = 15, units = "in", res = 500)
extractConsequenceSummary(
    exampleSwitchListAnalyzed_fet_sex,
    consequencesToAnalyze='all',
    plotGenes = FALSE,           # enables analysis of genes (instead of isoforms)
    asFractionTotal = FALSE)      # enables analysis of fraction of significant features
dev.off()

png("./Isoform_Consequence_Enrichment_fet_sex.png", height = 9, width = 15, units = "in", res = 500)
extractConsequenceEnrichment(
    exampleSwitchListAnalyzed_fet_sex,
    consequencesToAnalyze='all',
    analysisOppositeConsequence = TRUE,
    localTheme = theme_bw(base_size = 14), # Increase font size in vignette
    returnResult = FALSE # if TRUE returns a data.frame with the summary statistics
)
dev.off() #NO CONSEQUENCE ENRICHMENT FOR MY DATASET

#Overview Volcano Plot

## Each dot is an Isoform that has been switched (and sig and non-sig is indicated by low dIF)

png("./Isoform_Switch_Vol_all_genes_fet_sex.png", height = 9, width = 15, units = "in", res = 500)
ggplot(data=exampleSwitchListAnalyzed_fet_sex$isoformFeatures, aes(x=dIF, y=-log10(isoform_switch_q_value))) +
     geom_point(
        aes( color=abs(dIF) > 0.1 & isoform_switch_q_value < 0.05 ), # default cutoff
        size=1
    ) +
    geom_hline(yintercept = -log10(0.05), linetype='dashed') + # default cutoff
    geom_vline(xintercept = c(-0.1, 0.1), linetype='dashed') + # default cutoff for effect size
    facet_wrap( ~ condition_2) +
    #facet_grid(condition_1 ~ condition_2) + # alternative to facet_wrap if you have overlapping conditions
    scale_color_manual('Signficant\nIsoform Switch', values = c('black','red')) +
    labs(x='dIF', y='-Log10 ( Isoform Switch Q Value )') +
    theme_bw()
dev.off() 


#Switch vs Gene volcano plot (if gene expression changes and isoform switches are mutually exclusive)

png("./Isoform_Switch_vs_Gene_Vol_all_genes_fet_sex.png", height = 9, width = 15, units = "in", res = 500)
ggplot(data=exampleSwitchListAnalyzed_fet_sex$isoformFeatures, aes(x=gene_log2_fold_change, y=dIF)) +
    geom_point(
        aes( color=abs(dIF) > 0.1 & isoform_switch_q_value < 0.05 ), # default cutoff
        size=1
    ) + 
    facet_wrap(~ condition_2) +
    #facet_grid(condition_1 ~ condition_2) + # alternative to facet_wrap if you have overlapping conditions
    geom_hline(yintercept = 0, linetype='dashed') +
    geom_vline(xintercept = 0, linetype='dashed') +
    scale_color_manual('Signficant\nIsoform Switch', values = c('black','red')) +
    labs(x='Gene log2 fold change', y='dIF') +
    theme_bw()
dev.off() 







#Subset to Lipid Candidate Genes - Run Lipid Candidate Isoform Switch

lipid_ensembl_gencode <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/DESeq_rerun/lipid_ensembl_gencode.csv")

aSwitchList_Lip_fet_sex <- subsetSwitchAnalyzeRlist(
    switchAnalyzeRlist = aSwitchList_fet_sex,
    subset = aSwitchList_fet_sex$isoformFeatures$gene_name %in% lipid_ensembl_gencode$gene_symbol
)


#Identify Isoform Switches - Whole Genome

dir.create("./Isoform_DEXSeq_Lip_fet_sex")
aSwitchList_Lip_fet_sex <- preFilter(aSwitchList_Lip_fet_sex,
                isoCount = 10,
                min.Count.prop = 0.7,
                IFcutoff = 0.1,
                min.IF.prop = 0.5,
                alpha = 0.05,
                dIFcutoff = 0.1) 

exampleSwitchListAnalyzed_Lip_fet_sex <- isoformSwitchTestDEXSeq(
    switchAnalyzeRlist   = aSwitchList_Lip_fet_sex,
    reduceToSwitchingGenes = FALSE
)

#Summary of isoform switches
Isoform_Switch_Summary_Lip_fet_sex <- extractSwitchSummary(exampleSwitchListAnalyzed_Lip_fet_sex)
write.csv(Isoform_Switch_Summary_Lip_fet_sex, "Isoform_Switch_Summary_Lip_fet_sex.csv")

exampleSwitchListAnalyzed_Lip_fet_sex <- extractSequence(
    exampleSwitchListAnalyzed_Lip_fet_sex, 
    pathToOutput = "./Isoform_DEXSeq_Lip_fet_sex"
)

#In browser paste isoformSwitchAnalyzeR_isoform_nt.fasta file into https://cpc2.gao-lab.org/ to get information about coding potential (CPC2 - Coding Potential Calculator 2)

#In browser paste isoformSwitchAnalyzeR_isoform_AA.fasta file into https://www.ebi.ac.uk/jdispatcher/pfa/pfamscan to get information about protein families impacted

#In browser paste isoformSwitchAnalyzeR_isoform_AA.fasta file into https://iupred2a.elte.hu/ to predict switch to intrinsically disordered proteins

#In browser paste isoformSwitchAnalyzeR_isoform_AA.fasta file into https://services.healthtech.dtu.dk/services/SignalP-6.0/ to predict signal peptides and their cleavage sites


dir.create("./Isoform_P2_Lip_fet_sex")

##Copy all those results .txt files into "./Isoform_P2" directory

#Analysis of Alternative Splicing and functional consequences

exampleSwitchListAnalyzed_Lip_fet_sex <- isoformSwitchAnalysisPart2(
  switchAnalyzeRlist        = exampleSwitchListAnalyzed_Lip_fet_sex, 
  n                         = 50,    # if plotting was enabled, it would only output the top 10 switches
  removeNoncodinORFs        = FALSE,
  pathToCPC2resultFile      = "./Isoform_DEXSeq_Lip_fet_sex/cpc2_result.txt", #adds column identifying if transcript is protein- or non-coding
  pathToPFAMresultFile      = "./Isoform_DEXSeq_Lip_fet_sex/pfam_results.txt", #identifies functional and structural domains located on exons
  pathToIUPred2AresultFile  = "./Isoform_DEXSeq_Lip_fet_sex/iupred2a_result.txt", #marks disordered regions (involved in cell signaling)
  pathToSignalPresultFile   = "./Isoform_DEXSeq_Lip_fet_sex/signalP_results.txt", 
  outputPlots               = TRUE)

exampleSwitchListAnalyzed_Lip_fet_sex <- analyzeAlternativeSplicing(
  switchAnalyzeRlist = exampleSwitchListAnalyzed_Lip_fet_sex,
  quiet = FALSE
)
# write.csv(exampleSwitchList$isoformFeatures,"PE_isoform_features_Lip_fet_sex.csv", row.names = FALSE)

#Only females had sig_diff isoforms and Switches
 female_SwitchList <- subsetSwitchAnalyzeRlist(
  switchAnalyzeRlist = exampleSwitchListAnalyzed_Lip_fet_sex,
  subset = exampleSwitchListAnalyzed_Lip_fet_sex$isoformFeatures$condition_1 == "Control_F" &
           exampleSwitchListAnalyzed_Lip_fet_sex$isoformFeatures$condition_2 == "PE_F"
)

# female_SwitchList <- analyzeAlternativeSplicing(
#   switchAnalyzeRlist = female_SwitchList,
#   quiet = FALSE
# )

saveRDS(exampleSwitchListAnalyzed_Lip_fet_sex, file = "./Isoform_DEXSeq_Lip/exampleSwitchListAnalyzed_Lip_fet_sex.rds")

dir.create("./Lip_F_Isoform_Tables")
setwd("./Lip_F_Isoform_Tables")

Isoform_Table(female_SwitchList, "Lip_F", RNA_Seq_Lip_isoform_female, RNA_Lip_DESeq_results_autosomes_F, rMATS_results_WG)


#Plots

png("./Isoform_Splicing_Summary_Lip_fet_sex.png", height = 9, width = 15, units = "in", res = 500)
extractSplicingSummary(
    exampleSwitchListAnalyzed_Lip_fet_sex,
    asFractionTotal = FALSE,
    plotGenes=FALSE)
dev.off()

png("./Isoform_Splicing_Enrichment_Lip_fet_sex.png", height = 9, width = 15, units = "in", res = 500)
extractSplicingEnrichment(
    exampleSwitchListAnalyzed_Lip_fet_sex,
    returnResult = TRUE # if TRUE returns a data.frame with the summary statistics
)
dev.off() # NO SPLICING ENRICHMENT FOR MY DATASET



png("./Isoform_Switching_Summary_Lip_fet_sex.png", height = 9, width = 15, units = "in", res = 500)
extractTopSwitches(exampleSwitchListAnalyzed_Lip_fet_sex, filterForConsequences = TRUE)
dev.off()


#Consequence of Isoform Switches
png("./Isoform_Consequence_Summary_Lip_fet_sex.png", height = 9, width = 15, units = "in", res = 500)
extractConsequenceSummary(
    exampleSwitchListAnalyzed_Lip_fet_sex,
    consequencesToAnalyze='all',
    plotGenes = FALSE,           # enables analysis of genes (instead of isoforms)
    asFractionTotal = FALSE)      # enables analysis of fraction of significant features
dev.off()

png("./Isoform_Consequence_Enrichment_Lip_fet_sex.png", height = 9, width = 15, units = "in", res = 500)
extractConsequenceEnrichment(
    exampleSwitchListAnalyzed_Lip_fet_sex,
    consequencesToAnalyze='all',
    analysisOppositeConsequence = TRUE,
    localTheme = theme_bw(base_size = 14), # Increase font size in vignette
    returnResult = FALSE # if TRUE returns a data.frame with the summary statistics
)
dev.off() #NO CONSEQUENCE ENRICHMENT FOR MY DATASET

#Overview Volcano Plot

## Each dot is an Isoform that has been switched (and sig and non-sig is indicated by low dIF)

png("./Isoform_Switch_Vol_Lip_fet_sex.png", height = 9, width = 15, units = "in", res = 500)
ggplot(data=exampleSwitchListAnalyzed_Lip_fet_sex$isoformFeatures, aes(x=dIF, y=-log10(isoform_switch_q_value))) +
     geom_point(
        aes( color=abs(dIF) > 0.1 & isoform_switch_q_value < 0.05 ), # default cutoff
        size=1
    ) +
    geom_hline(yintercept = -log10(0.05), linetype='dashed') + # default cutoff
    geom_vline(xintercept = c(-0.1, 0.1), linetype='dashed') + # default cutoff for effect size
    facet_wrap( ~ condition_2) +
    #facet_grid(condition_1 ~ condition_2) + # alternative to facet_wrap if you have overlapping conditions
    scale_color_manual('Signficant\nIsoform Switch', values = c('black','red')) +
    labs(x='dIF', y='-Log10 ( Isoform Switch Q Value )') +
    theme_bw()
dev.off() 


#Switch vs Gene volcano plot (if gene expression changes and isoform switches are mutually exclusive)

png("./Isoform_Switch_vs_Gene_Vol_Lip_fet_sex.png", height = 9, width = 15, units = "in", res = 500)
ggplot(data=exampleSwitchListAnalyzed_Lip_fet_sex$isoformFeatures, aes(x=gene_log2_fold_change, y=dIF)) +
    geom_point(
        aes( color=abs(dIF) > 0.1 & isoform_switch_q_value < 0.05 ), # default cutoff
        size=1
    ) + 
    facet_wrap(~ condition_2) +
    #facet_grid(condition_1 ~ condition_2) + # alternative to facet_wrap if you have overlapping conditions
    geom_hline(yintercept = 0, linetype='dashed') +
    geom_vline(xintercept = 0, linetype='dashed') +
    scale_color_manual('Signficant\nIsoform Switch', values = c('black','red')) +
    labs(x='Gene log2 fold change', y='dIF') +
    theme_bw()
dev.off() 


rMATS_results_WG <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2026_Alternative_Splicing/IsoformSwitchAnalyzeR/rMATS_results_WG.csv")

#Read Differential Isoform results
Swish_isoform_combined <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2026_Alternative_Splicing/Swish/Swish_isoform_combined.csv")
Swish_isoform_female <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2026_Alternative_Splicing/Swish/Swish_isoform_female.csv")
Swish_isoform_male <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2026_Alternative_Splicing/Swish/Swish_isoform_male.csv")
Swish_Lip_isoform_combined <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2026_Alternative_Splicing/Swish/Swish_Lip_isoform_combined.csv")
Swish_Lip_isoform_female <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2026_Alternative_Splicing/Swish/Swish_Lip_isoform_female.csv")
Swish_Lip_isoform_male <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2026_Alternative_Splicing/Swish/Swish_Lip_isoform_male.csv")

#Read DESeq DEG results
RNA_WG_DESeq_results_autosomes_combined_sex <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_WG_rerun/RNA_WG_DESeq_results_autosomes_combined_sex.csv")
RNA_Lip_DESeq_results_autosomes_combined_sex <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/DESeq_rerun/RNA_Lip_DESeq_results_autosomes_combined_sex.csv")
RNA_Lip_sex_stratified <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/DESeq_rerun/RNA_Lip_DESeq_results_autosomes_IT_all.csv")
RNA_Lip_DESeq_results_autosomes_F <- RNA_Lip_sex_stratified[RNA_Lip_sex_stratified$comparison == "F_PEvsF_Cont", ]
RNA_Lip_DESeq_results_autosomes_M <- RNA_Lip_sex_stratified[RNA_Lip_sex_stratified$comparison == "M_PEvsM_Cont", ]
RNA_WG_sex_stratified <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/DESeq_rerun/RNA_WG_DESeq_results_autosomes_IT_all.csv")
RNA_WG_DESeq_results_autosomes_F <- RNA_WG_sex_stratified[RNA_WG_sex_stratified$comparison == "F_PEvsF_Cont", ]
RNA_WG_DESeq_results_autosomes_M <- RNA_WG_sex_stratified[RNA_WG_sex_stratified$comparison == "M_PEvsM_Cont", ]


#Make formation of tables into a function


Isoform_Table <- function (SwitchList, file_prefix, Swish_iso_name, DESeq2_name, rMATS_name){
isoformFeatures_filename <- paste0(file_prefix, "_", deparse(substitute(SwitchList)),".csv")
write.csv(SwitchList$isoformFeatures, isoformFeatures_filename)

condense_isoform <- SwitchList$isoformFeatures[SwitchList$isoformFeatures$isoform_switch_q_value <= 0.05 & (SwitchList$isoformFeatures$dIF <= -0.1 | SwitchList$isoformFeatures$dIF >= 0.1),]
sig_genes_isoform <- condense_isoform$gene_name
all_condensed_isoform <- SwitchList$isoformFeatures[SwitchList$isoformFeatures$gene_name %in% sig_genes_isoform,]

#Alternative Splicing
Isoform_Switch <- SwitchList$AlternativeSplicingAnalysis
Isoform_Switch$alt_splicing_type <- apply(Isoform_Switch[, c("ES", "MEE", "MES", "IR", "A5", "A3", "ATSS", "ATTS")], 1, function(row) {
  paste(names(row)[row == 1 & !is.na(row)], collapse = ";")
})

iso_alt <- merge(all_condensed_isoform, Isoform_Switch[,c("isoform_id","alt_splicing_type")], by = "isoform_id")

filename_isoform <- paste0(file_prefix, "_PE_Isoform_Feature.csv")
write.csv(iso_alt, filename_isoform, row.names = FALSE)

#My final table of switches
top_switches <- extractTopSwitches(SwitchList, filterForConsequences = FALSE, n = Inf)

up   <- subset(iso_alt, dIF > 0)
down <- subset(iso_alt, dIF < 0)

top_switches$isoformUpregulated <- tapply(up$isoform_id, up$gene_id, paste, collapse = "; ")[top_switches$gene_id]
top_switches$dIF_upregulated <- tapply(up$dIF,up$gene_id,paste, collapse = "; ")[top_switches$gene_id]
top_switches$alt_splic_upregulated <- tapply(up$alt_splicing_type,up$gene_id,paste, collapse = "; ")[top_switches$gene_id]

top_switches$isoformDownregulated <- tapply(down$isoform_id, down$gene_id, paste, collapse = "; ")[top_switches$gene_id]
top_switches$dIF_downregulated <- tapply(down$dIF, down$gene_id, paste, collapse = "; ")[top_switches$gene_id]
top_switches$alt_splic_down <- tapply(down$alt_splicing_type, down$gene_id, paste, collapse = "; ")[top_switches$gene_id]

switch_consequences <- SwitchList$switchConsequence
switch_consequences <- subset(SwitchList$switchConsequence, isoformsDifferent == TRUE & !is.na(switchConsequence))
switch_consequences_agg <- aggregate(switchConsequence ~ gene_name, data = switch_consequences, paste, collapse = "; ")
top_switches <- merge(top_switches, switch_consequences_agg, by = "gene_name", all.x = TRUE)
colnames(top_switches)[colnames(top_switches) == "switchConsequence"] <- "switchConsequences_type"

filename_top_switch<- paste0(file_prefix, "_PE_Isoform_Top_Switches.csv")
write.csv(top_switches, filename_top_switch)


#Make Differential Isoform (Swish) and Isoform Switch (IsoformSwitchAnalyzeR) merged table WG
Swish_top_switches <- top_switches
swish_lookup <- setNames(as.character(round(Swish_iso_name$Log2_FoldChg, 3)), Swish_iso_name$Transcript_ID)
Swish_top_switches$isoformUpregulated_Sig_Swish_Log2FC <- str_replace_all(Swish_top_switches$isoformUpregulated, swish_lookup)
Swish_top_switches$isoformDownregulated_Sig_Swish_Log2FC <- str_replace_all(Swish_top_switches$isoformDownregulated, swish_lookup)
# Swish_top_switches$Gene_ID.value <- Swish_iso_name$Gene_ID.value[match(Swish_top_switches$isoformDownregulated, Swish_iso_name$Transcript_ID)]
Swish_top_switches$Gene_ID.value <- gtf_data$gene_id[match(sub(";.*", "", Swish_top_switches$isoformUpregulated), gtf_data$transcript_id)]
Swish_top_switches$GeneID <- sub("\\..*", "", Swish_top_switches$Gene_ID.value)

#Make Differential Isoform (Swish) and Isoform Switch (IsoformSwitchAnalyzeR) + Swish DGE + rMATS merged table WG
DESeq2_name$Row.names <- sub("\\..*", "", DESeq2_name$Row.names)
DTE_Switch_DGE <- merge(Swish_top_switches, DESeq2_name, by.x = "GeneID", by.y = "Row.names", all.x = TRUE)
rMATS_name$GeneID <- sub("\\..*", "", rMATS_name$GeneID)
DTE_Switch_DGE_rMATS <- merge(DTE_Switch_DGE, rMATS_name, by = "GeneID", all.x = TRUE)

filename_DTE_Switch_DGE_rMATS <- paste0(file_prefix, "_Swish_IsoformSwitch.csv")
write.csv(DTE_Switch_DGE_rMATS, filename_DTE_Switch_DGE_rMATS) #0


#Make Comparisons only between Swish, rMATS and DESeq2
Swish_iso_name$gene_name <- gtf_data$gene_name[match(Swish_iso_name$Transcript_ID, gtf_data$transcript_id)]
Swish_iso_name$GeneID <- sub("\\..*", "", Swish_iso_name$Gene_ID.value)
Swish_DEG <-  merge(Swish_iso_name, DESeq2_name, by.x = "GeneID", by.y = "Row.names", all = TRUE)
Swish_rMATS_DEG <- merge(Swish_DEG, rMATS_name, by = "GeneID", all = TRUE)

filename_Swish_rMATS_DEG <- paste0(file_prefix, "_Swish_rMATS_DEG.csv")
write.csv(Swish_rMATS_DEG, filename_Swish_rMATS_DEG) 

#Identify percent of DEGs that were due to a single versus multiple differential transcript changes
transcript_counts <- Swish_rMATS_DEG %>% group_by(GeneID) %>% summarise(transcript_per_gene = n_distinct(Transcript_ID), .groups = "drop")
transcript_counts <- as.data.frame(transcript_counts)
Swish_rMATS_DEG_transcript <- Swish_rMATS_DEG %>% left_join(transcript_counts, by = "GeneID")

# head(Swish_rMATS_DEG_transcript)

sig_diff_ts <- subset(Swish_rMATS_DEG_transcript, Q_Value < 0.05)
sig_downregulated_ts <- subset(sig_diff_ts, Log2_FoldChg < -1)
sig_upregulated_ts <- subset(sig_diff_ts, Log2_FoldChg > 1)

downregulated_transcript_counts <- sig_downregulated_ts %>% group_by(GeneID) %>% summarise(sig_downregulated_transcripts = n_distinct(Transcript_ID), .groups = "drop_last")
downregulated_transcript_counts <- as.data.frame(downregulated_transcript_counts)
Swish_rMATS_DEG_transcript_down <- Swish_rMATS_DEG_transcript %>% left_join(downregulated_transcript_counts, by = "GeneID")
Swish_rMATS_DEG_transcript_down$sig_downregulated_transcripts[is.na(Swish_rMATS_DEG_transcript_down$sig_downregulated_transcripts)] <- 0

upregulated_transcript_counts <- sig_upregulated_ts %>% group_by(GeneID) %>% summarise(sig_upregulated_transcripts = n_distinct(Transcript_ID), .groups = "drop")
upregulated_transcript_counts <- as.data.frame(upregulated_transcript_counts)
Swish_rMATS_DEG_transcript_all <- Swish_rMATS_DEG_transcript_down %>% left_join(upregulated_transcript_counts, by = "GeneID")
Swish_rMATS_DEG_transcript_all$sig_upregulated_transcripts[is.na(Swish_rMATS_DEG_transcript_all$sig_upregulated_transcripts)] <- 0

Swish_rMATS_DEG_transcript_all$total_diff_transcripts <- Swish_rMATS_DEG_transcript_all$sig_downregulated_transcripts + Swish_rMATS_DEG_transcript_all$sig_upregulated_transcripts
Swish_rMATS_DEG_transcript_all$percent_diff_transcripts <- (Swish_rMATS_DEG_transcript_all$total_diff_transcripts)/(Swish_rMATS_DEG_transcript_all$transcript_per_gene)

#Only looking at sig diff expressed genes (DESeq2)
Swish_rMATS_sigDEG_transcript_all <- subset(Swish_rMATS_DEG_transcript_all, padj < 0.05)
Swish_rMATS_sigDEG_transcript_all <- subset(Swish_rMATS_sigDEG_transcript_all, (log2FoldChange < -1 | log2FoldChange > 1))

Swish_rMATS_sigDEG_transcript_all$DTE_group <- "several_isoforms_different"
all_isoforms_different <- subset(Swish_rMATS_sigDEG_transcript_all, total_diff_transcripts == transcript_per_gene)
Swish_rMATS_sigDEG_transcript_all$DTE_group[Swish_rMATS_sigDEG_transcript_all$total_diff_transcripts == 0] <- "no_isoforms_different"
Swish_rMATS_sigDEG_transcript_all$DTE_group[Swish_rMATS_sigDEG_transcript_all$total_diff_transcripts == 1] <- "one_isoform_different"
Swish_rMATS_sigDEG_transcript_all$DTE_group[Swish_rMATS_sigDEG_transcript_all$total_diff_transcripts == Swish_rMATS_sigDEG_transcript_all$transcript_per_gene] <- "all_isoforms_different" 

filename_transcript_DEG <- paste0(file_prefix, "_transcript_DEG.csv")
write.csv(Swish_rMATS_sigDEG_transcript_all, filename_transcript_DEG) 

#Overview Volcano Plot

# Each dot is an Isoform that has been switched (and sig and non-sig is indicated by low dIF)

vol_plot_filename <- paste0("./", file_prefix, "_Isoform_Switch_Vol.png")

png(vol_plot_filename, height = 9, width = 15, units = "in", res = 500)
iso_vol <- ggplot(data=SwitchList$isoformFeatures, aes(x=dIF, y=-log10(isoform_switch_q_value))) +
     geom_point(
        aes( color=abs(dIF) > 0.1 & isoform_switch_q_value < 0.05 ), # default cutoff
        size=1
    ) +
    geom_hline(yintercept = -log10(0.05), linetype='dashed') + # default cutoff
    geom_vline(xintercept = c(-0.1, 0.1), linetype='dashed') + # default cutoff for effect size
    facet_wrap( ~ condition_2) +
    #facet_grid(condition_1 ~ condition_2) + # alternative to facet_wrap if you have overlapping conditions
    scale_color_manual('Signficant\nIsoform Switch', values = c('black','red')) +
    labs(x='dIF', y='-Log10 ( Isoform Switch Q Value )') +
    theme_bw()

print(iso_vol)
dev.off() 


#Switch vs Gene volcano plot (if gene expression changes and isoform switches are mutually exclusive)

switch_v_gene_filename <- paste0("./", file_prefix, "_Isoform_Switch_vs_Gene_Vol.png")
png(switch_v_gene_filename, height = 9, width = 15, units = "in", res = 500)
gene_vol <- ggplot(data=SwitchList$isoformFeatures, aes(x=gene_log2_fold_change, y=dIF)) +
    geom_point(
        aes( color=abs(dIF) > 0.1 & isoform_switch_q_value < 0.05 ), # default cutoff
        size=1
    ) + 
    facet_wrap(~ condition_2) +
    #facet_grid(condition_1 ~ condition_2) + # alternative to facet_wrap if you have overlapping conditions
    geom_hline(yintercept = 0, linetype='dashed') +
    geom_vline(xintercept = 0, linetype='dashed') +
    scale_color_manual('Signficant\nIsoform Switch', values = c('black','red')) +
    labs(x='Gene log2 fold change', y='dIF') +
    theme_bw()

print(gene_vol)
dev.off() 

}



dir.create("./WG_F_Isoform_Tables")
setwd("./WG_F_Isoform_Tables")
Isoform_Table(WG_female_SwitchList, "WG_F", Swish_isoform_female, RNA_WG_DESeq_results_autosomes_F, rMATS_results_WG)


setwd("../")
dir.create("./WG_M_Isoform_Tables")
setwd("./WG_M_Isoform_Tables")
Isoform_Table(WG_male_SwitchList, "WG_M", Swish_isoform_male, RNA_WG_DESeq_results_autosomes_M, rMATS_results_WG)


dir.create("./Lip_F_Isoform_Tables")
setwd("./Lip_F_Isoform_Tables")

Isoform_Table(female_SwitchList, "Lip_F", Swish_Lip_isoform_female, RNA_Lip_DESeq_results_autosomes_F, rMATS_results_WG)


#Making tables for Lip_M
Isoform_Table_Lip_M("Lip_M", Swish_Lip_isoform_male, RNA_Lip_DESeq_results_autosomes_M, rMATS_results_WG)

Isoform_Table_Lip_M <- function (file_prefix, Swish_iso_name, DESeq2_name, rMATS_name){

#Make Comparisons only between Swish, rMATS and DESeq2
Swish_iso_name$gene_name <- gtf_data$gene_name[match(Swish_iso_name$Transcript_ID, gtf_data$transcript_id)]
Swish_iso_name$GeneID <- sub("\\..*", "", Swish_iso_name$Gene_ID.value)
DESeq2_name$gene_id <- sub("\\..*", "", DESeq2_name$gene_id)
rMATS_name$GeneID <- sub("\\..*", "", rMATS_name$GeneID)
Swish_DEG <-  merge(Swish_iso_name, DESeq2_name, by.x = "GeneID", by.y = "gene_id", all = TRUE)
Swish_rMATS_DEG <- merge(Swish_DEG, rMATS_name, by = "GeneID", all = TRUE)

filename_Swish_rMATS_DEG <- paste0(file_prefix, "_Swish_rMATS_DEG.csv")
write.csv(Swish_rMATS_DEG, filename_Swish_rMATS_DEG) 

#Identify percent of DEGs that were due to a single versus multiple differential transcript changes
transcript_counts <- Swish_rMATS_DEG %>% group_by(GeneID) %>% summarise(transcript_per_gene = n_distinct(Transcript_ID), .groups = "drop")
transcript_counts <- as.data.frame(transcript_counts)
Swish_rMATS_DEG_transcript <- Swish_rMATS_DEG %>% left_join(transcript_counts, by = "GeneID")

# head(Swish_rMATS_DEG_transcript)

sig_diff_ts <- subset(Swish_rMATS_DEG_transcript, Q_Value < 0.05)
sig_downregulated_ts <- subset(sig_diff_ts, Log2_FoldChg < -1)
sig_upregulated_ts <- subset(sig_diff_ts, Log2_FoldChg > 1)

downregulated_transcript_counts <- sig_downregulated_ts %>% group_by(GeneID) %>% summarise(sig_downregulated_transcripts = n_distinct(Transcript_ID), .groups = "drop_last")
downregulated_transcript_counts <- as.data.frame(downregulated_transcript_counts)
Swish_rMATS_DEG_transcript_down <- Swish_rMATS_DEG_transcript %>% left_join(downregulated_transcript_counts, by = "GeneID")
Swish_rMATS_DEG_transcript_down$sig_downregulated_transcripts[is.na(Swish_rMATS_DEG_transcript_down$sig_downregulated_transcripts)] <- 0

upregulated_transcript_counts <- sig_upregulated_ts %>% group_by(GeneID) %>% summarise(sig_upregulated_transcripts = n_distinct(Transcript_ID), .groups = "drop")
upregulated_transcript_counts <- as.data.frame(upregulated_transcript_counts)
Swish_rMATS_DEG_transcript_all <- Swish_rMATS_DEG_transcript_down %>% left_join(upregulated_transcript_counts, by = "GeneID")
Swish_rMATS_DEG_transcript_all$sig_upregulated_transcripts[is.na(Swish_rMATS_DEG_transcript_all$sig_upregulated_transcripts)] <- 0

Swish_rMATS_DEG_transcript_all$total_diff_transcripts <- Swish_rMATS_DEG_transcript_all$sig_downregulated_transcripts + Swish_rMATS_DEG_transcript_all$sig_upregulated_transcripts
Swish_rMATS_DEG_transcript_all$percent_diff_transcripts <- (Swish_rMATS_DEG_transcript_all$total_diff_transcripts)/(Swish_rMATS_DEG_transcript_all$transcript_per_gene)

#Only looking at sig diff expressed genes (DESeq2)
Swish_rMATS_sigDEG_transcript_all <- subset(Swish_rMATS_DEG_transcript_all, padj < 0.05)
Swish_rMATS_sigDEG_transcript_all <- subset(Swish_rMATS_sigDEG_transcript_all, (log2FoldChange < -1 | log2FoldChange > 1))

Swish_rMATS_sigDEG_transcript_all$DTE_group <- "several_isoforms_different"
all_isoforms_different <- subset(Swish_rMATS_sigDEG_transcript_all, total_diff_transcripts == transcript_per_gene)
Swish_rMATS_sigDEG_transcript_all$DTE_group[Swish_rMATS_sigDEG_transcript_all$total_diff_transcripts == 0] <- "no_isoforms_different"
Swish_rMATS_sigDEG_transcript_all$DTE_group[Swish_rMATS_sigDEG_transcript_all$total_diff_transcripts == 1] <- "one_isoform_different"
Swish_rMATS_sigDEG_transcript_all$DTE_group[Swish_rMATS_sigDEG_transcript_all$total_diff_transcripts == Swish_rMATS_sigDEG_transcript_all$transcript_per_gene] <- "all_isoforms_different" 

filename_transcript_DEG <- paste0(file_prefix, "_transcript_DEG.csv")
write.csv(Swish_rMATS_sigDEG_transcript_all, filename_transcript_DEG) 

}

setwd("./Lip_M_Isoform_Tables")
Isoform_Table_Lip_M("Lip_M", Swish_Lip_isoform_male, RNA_Lip_DESeq_results_autosomes_M, rMATS_results_WG)


#Load CSVs for Stacked Bar Plot 
Lip_combined_stacked <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2026_Alternative_Splicing/IsoformSwitchAnalyzeR/Lip_combined_sex_Isoform_Tables/Lip_combined_sex_transcript_DEG.csv")
Lip_F_stacked <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2026_Alternative_Splicing/IsoformSwitchAnalyzeR/Lip_F_Isoform_Tables/Lip_F_transcript_DEG.csv")
Lip_M_stacked <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2026_Alternative_Splicing/IsoformSwitchAnalyzeR/Lip_F_Isoform_Tables/Lip_M_transcript_DEG.csv")
WG_combined_stacked <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2026_Alternative_Splicing/IsoformSwitchAnalyzeR/WG_combined_sex_Isoform_Tables/WG_combined_sex_transcript_DEG.csv")
WG_F_stacked <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2026_Alternative_Splicing/IsoformSwitchAnalyzeR/WG_F_Isoform_Tables/WG_F_transcript_DEG.csv")
WG_M_stacked <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2026_Alternative_Splicing/IsoformSwitchAnalyzeR/WG_M_Isoform_Tables/WG_M_transcript_DEG.csv")


library(dplyr)
library(ggplot2)
library(scales) # Needed for percentage formatting on the y-axis

# 1. Read the CSVs and add a 'Group' column to identify each dataset
Lip_combined_stacked <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2026_Alternative_Splicing/IsoformSwitchAnalyzeR/Lip_combined_sex_Isoform_Tables/Lip_combined_sex_transcript_DEG.csv") %>% 
  mutate(Group = "Lip_combined")

Lip_F_stacked <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2026_Alternative_Splicing/IsoformSwitchAnalyzeR/Lip_F_Isoform_Tables/Lip_F_transcript_DEG.csv") %>% 
  mutate(Group = "Lip_F")

Lip_M_stacked <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2026_Alternative_Splicing/IsoformSwitchAnalyzeR/Lip_F_Isoform_Tables/Lip_M_transcript_DEG.csv") %>% 
  mutate(Group = "Lip_M")

WG_combined_stacked <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2026_Alternative_Splicing/IsoformSwitchAnalyzeR/WG_combined_sex_Isoform_Tables/WG_combined_sex_transcript_DEG.csv") %>% 
  mutate(Group = "WG_combined")

WG_F_stacked <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2026_Alternative_Splicing/IsoformSwitchAnalyzeR/WG_F_Isoform_Tables/WG_F_transcript_DEG.csv") %>% 
  mutate(Group = "WG_F")

WG_M_stacked <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2026_Alternative_Splicing/IsoformSwitchAnalyzeR/WG_M_Isoform_Tables/WG_M_transcript_DEG.csv") %>% 
  mutate(Group = "WG_M")

# 2. Combine all dataframes into a single dataframe
all_groups_combined <- bind_rows(
  Lip_combined_stacked, 
  Lip_F_stacked, 
  Lip_M_stacked, 
  WG_combined_stacked, 
  WG_F_stacked, 
  WG_M_stacked
)

# 3. Reset X-axis to your original fixed order
all_groups_combined$Group <- factor(
  all_groups_combined$Group, 
  levels = c("Lip_combined", "Lip_F", "Lip_M", "WG_combined", "WG_F", "WG_M")
)

# 4. Explicitly set stacking order from bottom to top
all_groups_combined$DTE_group <- factor(
  all_groups_combined$DTE_group, 
  levels = c(
    "no_isoforms_different",
    "all_isoforms_different", 
    "several_isoforms_different", 
    "one_isoform_different"
  )
)

# Initialize the PNG graphics device
png(file = "./Percent_Transcripts_Driving_DEGs.png", height = 7, width = 10, units = "in", res = 300)

# 5. Create the 100% stacked bar plot
print(
  ggplot(all_groups_combined, aes(x = Group, fill = DTE_group)) +
    geom_bar(position = "fill") + 
    scale_y_continuous(labels = scales::percent) + 
    labs(
      title = "Proportion of DEGs Driven by Transcript Changes",
      x = "Dataset Group",
      y = "Percentage of DEGs",
      fill = "DTE Group"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
)

# Close the file device
dev.off()













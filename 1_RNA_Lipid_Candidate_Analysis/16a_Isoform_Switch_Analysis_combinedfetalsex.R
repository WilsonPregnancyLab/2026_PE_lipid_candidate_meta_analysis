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
colnames(metadata)[colnames(metadata) == "disease_group"] <- "condition"

columns <- c("sampleID", "GSE_number", "condition", "predicted_fetal_sex")
myDesign <- metadata[, colnames(metadata) %in% columns]
myDesign_ordered <- myDesign[, c("sampleID", "condition", "GSE_number", "predicted_fetal_sex")]
myDesign_ordered$sampleID   <- as.character(myDesign_ordered$sampleID)
myDesign_ordered$condition  <- as.character(myDesign_ordered$condition)
myDesign_ordered$GSE_number <- as.character(myDesign_ordered$GSE_number)
myDesign_ordered$predicted_fetal_sex <- as.character(myDesign_ordered$predicted_fetal_sex)

# all(file.exists(myDesign$files))

aSwitchList <- importRdata(
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

dir.create("./Isoform_DEXSeq")
aSwitchList <- preFilter(aSwitchList,
                isoCount = 10,
                min.Count.prop = 0.7,
                IFcutoff = 0.1,
                min.IF.prop = 0.5,
                alpha = 0.05,
                dIFcutoff = 0.1) 

exampleSwitchListAnalyzed <- isoformSwitchTestDEXSeq(
    switchAnalyzeRlist   = aSwitchList,
    reduceToSwitchingGenes = FALSE
)

#Summary of isoform switches
Isoform_Switch_Summary <- extractSwitchSummary(exampleSwitchListAnalyzed)
write.csv(Isoform_Switch_Summary, "Isoform_Switch_Summary_all_genes.csv")

exampleSwitchListAnalyzed <- extractSequence(
    exampleSwitchListAnalyzed, 
    pathToOutput = "./Isoform_DEXSeq"
)

#In browser paste isoformSwitchAnalyzeR_isoform_nt.fasta file into https://cpc2.gao-lab.org/ to get information about coding potential (CPC2 - Coding Potential Calculator 2)

#In browser paste isoformSwitchAnalyzeR_isoform_AA.fasta file into https://www.ebi.ac.uk/jdispatcher/pfa/pfamscan to get information about protein families impacted

#In browser paste isoformSwitchAnalyzeR_isoform_AA.fasta file into https://services.healthtech.dtu.dk/services/SignalP-6.0/ to predict signal peptides and their cleavage sites

#In browser paste isoformSwitchAnalyzeR_isoform_AA.fasta file into https://services.healthtech.dtu.dk/services/SignalP-6.0/ to predict signal peptides and their cleavage sites


dir.create("./Isoform_P2")

##Copy all those results .txt files into "./Isoform_P2" directory

#Analysis of Alternative Splicing and functional consequences

exampleSwitchListAnalyzed <- isoformSwitchAnalysisPart2(
  switchAnalyzeRlist        = exampleSwitchListAnalyzed, 
  n                         = 50,    # if plotting was enabled, it would only output the top 10 switches
  removeNoncodinORFs        = TRUE,
  pathToCPC2resultFile      = "./Isoform_P2/cpc2_result.txt", #adds column identifying if transcript is protein- or non-coding
  pathToPFAMresultFile      = "./Isoform_P2/pfam_results.txt", #identifies functional and structural domains located on exons
  pathToIUPred2AresultFile  = "./Isoform_P2/iupred2a_result.txt", 
  pathToSignalPresultFile   = "./Isoform_P2/signalP_results.txt", #marks disordered regions (involved in cell signaling)
  outputPlots               = TRUE)


names(exampleSwitchListAnalyzed)
 [1] "isoformFeatures"             "exons"
 [3] "conditions"                  "designMatrix"
 [5] "sourceId"                    "isoformCountMatrix"
 [7] "isoformRepExpression"        "runInfo"
 [9] "orfAnalysis"                 "isoformRepIF"
[11] "ntSequence"                  "isoformSwitchAnalysis"
[13] "aaSequence"                  "domainAnalysis"
[15] "idrAnalysis"                 "signalPeptideAnalysis"
[17] "AlternativeSplicingAnalysis" "switchConsequence"

saveRDS(exampleSwitchListAnalyzed, file = "./Isoform_DEXSeq/exampleSwitchListAnalyzed.rds")


#Make single table with all data
condense_isoform <- exampleSwitchListAnalyzed$isoformFeatures[exampleSwitchListAnalyzed$isoformFeatures$isoform_switch_q_value <= 0.05 & (exampleSwitchListAnalyzed$isoformFeatures$dIF <= -0.1 | exampleSwitchListAnalyzed$isoformFeatures$dIF >= 0.1),]
sig_genes_isoform <- condense_isoform$gene_name
all_condensed_isoform <- exampleSwitchListAnalyzed$isoformFeatures[exampleSwitchListAnalyzed$isoformFeatures$gene_name %in% sig_genes_isoform,]

#Alternative Splicing
Isoform_Switch <- exampleSwitchListAnalyzed$AlternativeSplicingAnalysis
Isoform_Switch$alt_splicing_type <- apply(Isoform_Switch[, c("ES", "MEE", "MES", "IR", "A5", "A3", "ATSS", "ATTS")], 1, function(row) {
  paste(names(row)[row == 1 & !is.na(row)], collapse = ";")
})

iso_alt <- merge(all_condensed_isoform, Isoform_Switch[,c("isoform_id","alt_splicing_type")], by = "isoform_id")

write.csv(iso_alt,"WG_combined_sex_PE_isoform_features.csv", row.names = FALSE)

# #My final table of switches
# top_switches <- extractTopSwitches(exampleSwitchListAnalyzed, filterForConsequences = FALSE)

# up   <- subset(iso_alt, dIF > 0)
# down <- subset(iso_alt, dIF < 0)

# top_switches$isoformUpregulated <- tapply(up$isoform_id, up$gene_id, paste, collapse = "; ")[top_switches$gene_id]
# top_switches$dIF_upregulated <- tapply(up$dIF,up$gene_id,paste, collapse = "; ")[top_switches$gene_id]
# top_switches$alt_splic_upregulated <- tapply(up$alt_splicing_type,up$gene_id,paste, collapse = "; ")[top_switches$gene_id]

# top_switches$isoformDownregulated <- tapply(down$isoform_id, down$gene_id, paste, collapse = "; ")[top_switches$gene_id]
# top_switches$dIF_downregulated <- tapply(down$dIF, down$gene_id, paste, collapse = "; ")[top_switches$gene_id]
# top_switches$alt_splic_down <- tapply(down$alt_splicing_type, down$gene_id, paste, collapse = "; ")[top_switches$gene_id]

# switch_consequences <- exampleSwitchListAnalyzed$switchConsequence
# switch_consequences <- subset(exampleSwitchListAnalyzed$switchConsequence, isoformsDifferent == TRUE & !is.na(switchConsequence))
# switch_consequences_agg <- aggregate(switchConsequence ~ gene_name, data = switch_consequences, paste, collapse = "; ")
# top_switches <- merge(top_switches, switch_consequences_agg, by = "gene_name", all.x = TRUE)
# colnames(top_switches)[colnames(top_switches) == "switchConsequence"] <- "switchConsequences_type"



# Swish_isoform_combined <- Swish_isoform_combined[Swish_isoform_combined$Q_Value <= 0.05, c("Transcript_ID", "Gene_ID.value", "Log2_FoldChg", "Q_Value")]

# #Make Differential Isoform (Swish) and Isoform Switch (IsoformSwitchAnalyzeR) merged table WG
# Swish_top_switches <- top_switches
# Swish_top_switches$isoformUpregulated_Sig_Swish_Log2FC <- Swish_isoform_combined$Log2_FoldChg[match(Swish_top_switches$isoformUpregulated, Swish_isoform_combined$Transcript_ID)]
# Swish_top_switches$isoformDownregulated_Sig_Swish_Log2FC <- Swish_isoform_combined$Log2_FoldChg[match(Swish_top_switches$isoformDownregulated, Swish_isoform_combined$Transcript_ID)]
# Swish_top_switches$Gene_ID.value <- Swish_isoform_combined$Gene_ID.value[match(sub(";.*", "", Swish_top_switches$isoformUpregulated), Swish_isoform_combined$Transcript_ID)]


# # write.csv(Swish_top_switches, filename_DTE_Switch)

# #Make Differential Isoform (Swish) and Isoform Switch (IsoformSwitchAnalyzeR) + DGE + rMATS merged table WG
# DTE_Switch_DGE <- merge(Swish_top_switches, RNA_WG_DESeq_results_autosomes_combined_sex, by.x = "Gene_ID.value", by.y = "Row.names", all.x = TRUE)
# # DTE_Switch_DGE$GeneID <- sub("\\..*", "", DTE_Switch_DGE$Gene_ID.value)
# DTE_Switch_DGE_rMATS <- merge(DTE_Switch_DGE, rMATS_results_WG, by.x = "Gene_ID.value", by.y = "GeneID", all.x = TRUE)

# filename_DTE_Switch_DGE_rMATS <- paste0(file_prefix, "_Swish_IsoformSwitch.csv")
# write.csv(DTE_Switch_DGE_rMATS, rMATS_results_WG) #0



#Plots

# png("./Isoform_Splicing_Summary.png", height = 9, width = 15, units = "in", res = 500)
# extractSplicingSummary(
#     exampleSwitchListAnalyzed,
#     asFractionTotal = FALSE,
#     plotGenes=FALSE)
# dev.off()

# png("./Isoform_Splicing_Enrichment.png", height = 9, width = 15, units = "in", res = 500)
# extractSplicingEnrichment(
#     exampleSwitchListAnalyzed,
#     returnResult = TRUE # if TRUE returns a data.frame with the summary statistics
# )
# dev.off() # NO SPLICING ENRICHMENT FOR MY DATASET


# png("./Isoform_Switching_Summary.png", height = 9, width = 15, units = "in", res = 500)
# extractTopSwitches(exampleSwitchListAnalyzed, filterForConsequences = TRUE)
# dev.off()

# #Consequence of Isoform Switches
# png("./Isoform_Consequence_Summary.png", height = 9, width = 15, units = "in", res = 500)
# extractConsequenceSummary(
#     exampleSwitchListAnalyzed,
#     consequencesToAnalyze='all',
#     plotGenes = FALSE,           # enables analysis of genes (instead of isoforms)
#     asFractionTotal = FALSE)      # enables analysis of fraction of significant features
# dev.off()

# png("./Isoform_Consequence_Enrichment.png", height = 9, width = 15, units = "in", res = 500)
# extractConsequenceEnrichment(
#     exampleSwitchListAnalyzed,
#     consequencesToAnalyze='all',
#     analysisOppositeConsequence = TRUE,
#     localTheme = theme_bw(base_size = 14), # Increase font size in vignette
#     returnResult = FALSE # if TRUE returns a data.frame with the summary statistics
# )
# dev.off() #NO CONSEQUENCE ENRICHMENT FOR MY DATASET

#Overview Volcano Plot

## Each dot is an Isoform that has been switched (and sig and non-sig is indicated by low dIF)

png("./Isoform_Switch_Vol_all_genes.png", height = 9, width = 15, units = "in", res = 500)
ggplot(data=exampleSwitchListAnalyzed$isoformFeatures, aes(x=dIF, y=-log10(isoform_switch_q_value))) +
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

png("./Isoform_Switch_vs_Gene_Vol_all_genes.png", height = 9, width = 15, units = "in", res = 500)
ggplot(data=exampleSwitchListAnalyzed$isoformFeatures, aes(x=gene_log2_fold_change, y=dIF)) +
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

aSwitchList_Lip <- subsetSwitchAnalyzeRlist(
    switchAnalyzeRlist = aSwitchList,
    subset = aSwitchList$isoformFeatures$gene_name %in% lipid_ensembl_gencode$gene_symbol
)


#Identify Isoform Switches - Whole Genome

dir.create("./Isoform_DEXSeq_Lip")
aSwitchList_Lip <- preFilter(aSwitchList_Lip,
                isoCount = 10,
                min.Count.prop = 0.7,
                IFcutoff = 0.1,
                min.IF.prop = 0.5,
                alpha = 0.05,
                dIFcutoff = 0.1) 

exampleSwitchListAnalyzed_Lip <- isoformSwitchTestDEXSeq(
    switchAnalyzeRlist   = aSwitchList_Lip,
    reduceToSwitchingGenes = FALSE
)


#Summary of isoform switches
Isoform_Switch_Summary_Lip <- extractSwitchSummary(exampleSwitchListAnalyzed_Lip)
write.csv(Isoform_Switch_Summary_Lip, "Isoform_Switch_Summary_Lip.csv")

exampleSwitchListAnalyzed_Lip <- extractSequence(
    exampleSwitchListAnalyzed_Lip, 
    pathToOutput = "./Isoform_DEXSeq_Lip"
)

#In browser paste isoformSwitchAnalyzeR_isoform_nt.fasta file into https://cpc2.gao-lab.org/ to get information about coding potential (CPC2 - Coding Potential Calculator 2)

#In browser paste isoformSwitchAnalyzeR_isoform_AA.fasta file into https://www.ebi.ac.uk/jdispatcher/pfa/pfamscan to get information about protein families impacted

#In browser paste isoformSwitchAnalyzeR_isoform_AA.fasta file into https://iupred2a.elte.hu/ to predict switch to intrinsically disordered proteins

#In browser paste isoformSwitchAnalyzeR_isoform_AA.fasta file into https://services.healthtech.dtu.dk/services/SignalP-6.0/ to predict signal peptides and their cleavage sites


dir.create("./Isoform_P2_Lip")

##Copy all those results .txt files into "./Isoform_P2" directory

#Analysis of Alternative Splicing and functional consequences

exampleSwitchListAnalyzed_Lip <- isoformSwitchAnalysisPart2(
  switchAnalyzeRlist        = exampleSwitchListAnalyzed_Lip, 
  n                         = 50,    # if plotting was enabled, it would only output the top 10 switches
  removeNoncodinORFs        = TRUE,
  pathToCPC2resultFile      = "./Isoform_DEXSeq_Lip/cpc2_result.txt", #adds column identifying if transcript is protein- or non-coding
  pathToPFAMresultFile      = "./Isoform_DEXSeq_Lip/pfam_results.txt", #identifies functional and structural domains located on exons
  pathToIUPred2AresultFile  = "./Isoform_DEXSeq_Lip/iupred2a_result.txt", #marks disordered regions (involved in cell signaling)
  pathToSignalPresultFile   = "./Isoform_DEXSeq_Lip/signalP_results.txt", 
  outputPlots               = TRUE)

write.csv(exampleSwitchListAnalyzed_Lip$isoformFeatures,"./Isoform_DEXSeq_Lip/PE_isoform_features_Lip.csv", row.names = FALSE)

saveRDS(exampleSwitchListAnalyzed_Lip, file = "./Isoform_DEXSeq_Lip/exampleSwitchListAnalyzed_Lip.rds")


#Plots
setwd("./Isoform_DEXSeq_Lip/")

png("./Isoform_Splicing_Summary_Lip.png", height = 9, width = 15, units = "in", res = 500)
extractSplicingSummary(
    exampleSwitchListAnalyzed_Lip,
    asFractionTotal = FALSE,
    plotGenes=FALSE)
dev.off()

png("./Isoform_Splicing_Enrichment_Lip.png", height = 9, width = 15, units = "in", res = 500)
extractSplicingEnrichment(
    exampleSwitchListAnalyzed_Lip,
    returnResult = TRUE # if TRUE returns a data.frame with the summary statistics
)
dev.off() # NO SPLICING ENRICHMENT FOR MY DATASET



png("./Isoform_Switching_Summary_Lip.png", height = 9, width = 15, units = "in", res = 500)
extractTopSwitches(exampleSwitchListAnalyzed_Lip, filterForConsequences = TRUE)
dev.off()


#Consequence of Isoform Switches
png("./Isoform_Consequence_Summary_Lip.png", height = 9, width = 15, units = "in", res = 500)
extractConsequenceSummary(
    exampleSwitchListAnalyzed_Lip,
    consequencesToAnalyze='all',
    plotGenes = FALSE,           # enables analysis of genes (instead of isoforms)
    asFractionTotal = FALSE)      # enables analysis of fraction of significant features
dev.off()

png("./Isoform_Consequence_Enrichment_Lip.png", height = 9, width = 15, units = "in", res = 500)
extractConsequenceEnrichment(
    exampleSwitchListAnalyzed_Lip,
    consequencesToAnalyze='all',
    analysisOppositeConsequence = TRUE,
    localTheme = theme_bw(base_size = 14), # Increase font size in vignette
    returnResult = FALSE # if TRUE returns a data.frame with the summary statistics
)
dev.off() #NO CONSEQUENCE ENRICHMENT FOR MY DATASET

#Overview Volcano Plot

## Each dot is an Isoform that has been switched (and sig and non-sig is indicated by low dIF)

png("./Isoform_Switch_Vol_Lip.png", height = 9, width = 15, units = "in", res = 500)
ggplot(data=exampleSwitchListAnalyzed_Lip$isoformFeatures, aes(x=dIF, y=-log10(isoform_switch_q_value))) +
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

png("./Isoform_Switch_vs_Gene_Vol_Lip.png", height = 9, width = 15, units = "in", res = 500)
ggplot(data=exampleSwitchListAnalyzed_Lip$isoformFeatures, aes(x=gene_log2_fold_change, y=dIF)) +
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



#Comparing Alternative Splicing to Isoform
library(dplyr)
library(purrr)

setwd("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2026_Alternative_Splicing")

#Read alternative splicing results and combine to one table

SC_combined_JC <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2026_Alternative_Splicing/alt_splic_analysis/output/SE.MATS.JC.csv") #Skipped Exon (positive is more skipping in disease)
SC_combined_subset <- SC_combined_JC[,c("GeneID", "geneSymbol", "chr", "exonStart_0base", "exonEnd", "FDR", "IncLevelDifference")]
SC_combined_subset <- SC_combined_subset[SC_combined_subset$FDR <= 0.05 & (SC_combined_subset$IncLevelDifference <= -0.1 | SC_combined_subset$IncLevelDifference >= 0.1),]
colnames(SC_combined_subset)[colnames(SC_combined_subset) == "FDR"] <- "SC_FDR"
colnames(SC_combined_subset)[colnames(SC_combined_subset) == "IncLevelDifference"] <- "SC_IncLevelDifference"
colnames(SC_combined_subset)[colnames(SC_combined_subset) == "exonStart_0base"] <- "SC_exonStart_0base"
colnames(SC_combined_subset)[colnames(SC_combined_subset) == "exonEnd"] <- "SC_exonEnd"

RI_combined_JC <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2026_Alternative_Splicing/alt_splic_analysis/output/RI.MATS.JC.csv") #Retained Intron (positive is intron is more present in disease)
RI_combined_subset <- RI_combined_JC[,c("GeneID", "riExonStart_0base", "riExonEnd", "FDR", "IncLevelDifference")]
RI_combined_subset <- RI_combined_subset[RI_combined_subset$FDR <= 0.05 & (RI_combined_subset$IncLevelDifference <= -0.1 | RI_combined_subset$IncLevelDifference >= 0.1),]
colnames(RI_combined_subset)[colnames(RI_combined_subset) == "FDR"] <- "RI_FDR"
colnames(RI_combined_subset)[colnames(RI_combined_subset) == "IncLevelDifference"] <- "RI_IncLevelDifference"

MXE_combined_JC <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2026_Alternative_Splicing/alt_splic_analysis/output/MXE.MATS.JC.csv") #Mutually Exclusive Exons (positive is exon 1 is more often in controls, exon 2 is more often in PE)
MXE_combined_subset <- MXE_combined_JC[,c("GeneID", "X1stExonStart_0base", "X1stExonEnd", "X2ndExonStart_0base", "X2ndExonEnd", "FDR", "IncLevelDifference")]
MXE_combined_subset <- MXE_combined_subset[MXE_combined_subset$FDR <= 0.05 & (MXE_combined_subset$IncLevelDifference <= -0.1 | MXE_combined_subset$IncLevelDifference >= 0.1),]
colnames(MXE_combined_subset)[colnames(MXE_combined_subset) == "FDR"] <- "MXE_FDR"
colnames(MXE_combined_subset)[colnames(MXE_combined_subset) == "IncLevelDifference"] <- "MXE_IncLevelDifference"
colnames(MXE_combined_subset)[colnames(MXE_combined_subset) == "X1stExonStart_0base"] <- "MXE_X1stExonStart_0base"
colnames(MXE_combined_subset)[colnames(MXE_combined_subset) == "X1stExonEnd"] <- "MXE_X1stExonEnd"
colnames(MXE_combined_subset)[colnames(MXE_combined_subset) == "X2ndExonStart_0base"] <- "MXE_X2ndExonStart_0base"
colnames(MXE_combined_subset)[colnames(MXE_combined_subset) == "X2ndExonEnd"] <- "MXE_X2ndExonEnd"

A5SS_combined_JC <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2026_Alternative_Splicing/alt_splic_analysis/output/A5SS.MATS.JC.csv") #Alternative 5' Splice Site (positive is 5' splice site is shaved off on PE)
A5SS_combined_subset <- A5SS_combined_JC[,c("GeneID", "longExonStart_0base", "longExonEnd", "FDR", "IncLevelDifference")]
A5SS_combined_subset <- A5SS_combined_subset[A5SS_combined_subset$FDR <= 0.05 & (A5SS_combined_subset$IncLevelDifference <= -0.1 | A5SS_combined_subset$IncLevelDifference >= 0.1),]
colnames(A5SS_combined_subset)[colnames(A5SS_combined_subset) == "FDR"] <- "A5SS_FDR"
colnames(A5SS_combined_subset)[colnames(A5SS_combined_subset) == "IncLevelDifference"] <- "A5SS_IncLevelDifference"
colnames(A5SS_combined_subset)[colnames(A5SS_combined_subset) == "longExonStart_0base"] <- "A5SS_longExonStart_0base"
colnames(A5SS_combined_subset)[colnames(A5SS_combined_subset) == "longExonEnd"] <- "A5SS_longExonEnd"

A3SS_combined_JC <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2026_Alternative_Splicing/alt_splic_analysis/output/A3SS.MATS.JC.csv") #Alternative 3' Splicing Site (positive is 3' splice site is shaved off on PE)
A3SS_combined_subset <- A3SS_combined_JC[,c("GeneID", "longExonStart_0base", "longExonEnd", "FDR", "IncLevelDifference")]
A3SS_combined_subset <- A3SS_combined_subset[A3SS_combined_subset$FDR <= 0.05 & (A3SS_combined_subset$IncLevelDifference <= -0.1 | A3SS_combined_subset$IncLevelDifference >= 0.1),]
colnames(A3SS_combined_subset)[colnames(A3SS_combined_subset) == "FDR"] <- "A3SS_FDR"
colnames(A3SS_combined_subset)[colnames(A3SS_combined_subset) == "IncLevelDifference"] <- "A3SS_IncLevelDifference"
colnames(A3SS_combined_subset)[colnames(A3SS_combined_subset) == "longExonStart_0base"] <- "A3SS_longExonStart_0base"
colnames(A3SS_combined_subset)[colnames(A3SS_combined_subset) == "longExonEnd"] <- "A3SS_longExonEnd"

rMATS_list <- list(SC_combined_subset, RI_combined_subset, MXE_combined_subset, A5SS_combined_subset, A3SS_combined_subset)
rMATS_results_WG <- rMATS_list %>% reduce(full_join, by = "GeneID")
write.csv(rMATS_results_WG, "rMATS_results_WG.csv")

rMATS_results_WG <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2026_Alternative_Splicing/IsoformSwitchAnalyzeR/rMATS_results_WG.csv")

#Read Differential Isoform results
Swish_isoform_combined <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2026_Alternative_Splicing/Swish/Swish_isoform_combined.csv")
Swish_isoform_female <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2026_Alternative_Splicing/Swish/Swish_isoform_female.csv")
Swish_isoform_male <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2026_Alternative_Splicing/Swish/Swish_isoform_male.csv")
Swish_Lip_isoform_combined <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2026_Alternative_Splicing/Swish/Swish_Lip_isoform_combined.csv")
Swish_Lip_isoform_female <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2026_Alternative_Splicing/Swish/Swish_Lip_isoform_female.csv")
Swish_Lip_isoform_male <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2026_Alternative_Splicing/Swish/Swish_Lip_isoform_male.csv")

#Read Swish DEG results
Swish_DEG_combined <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2026_Alternative_Splicing/Swish/Swish_DEG_combined.csv")
Swish_DEG_F <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2026_Alternative_Splicing/Swish/Swish_DEG_F.csv")
Swish_DEG_M <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2026_Alternative_Splicing/Swish/Swish_DEG_M.csv")
Swish_DEG_Lip_combined <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2026_Alternative_Splicing/Swish/Swish_DEG_Lip_combined.csv")
Swish_DEG_Lip_F <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2026_Alternative_Splicing/Swish/Swish_DEG_Lip_F.csv")
Swish_DEG_Lip_M <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2026_Alternative_Splicing/Swish/Swish_DEG_Lip_M.csv")


#Read Isoform Switch results
#Isoform_Switch_WG <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2026_Alternative_Splicing/IsoformSwitchAnalyzeR/WG_combined_sex_PE_Isoform_Gene_Consequences.csv")

#Read DESeq DEG results
RNA_WG_DESeq_results_autosomes_combined_sex <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_WG_rerun/RNA_WG_DESeq_results_autosomes_combined_sex.csv")
RNA_Lip_DESeq_results_autosomes_combined_sex <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/DESeq_rerun/RNA_Lip_DESeq_results_autosomes_combined_sex.csv")
RNA_Lip_sex_stratified <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/DESeq_rerun/RNA_Lip_DESeq_results_autosomes_IT_all.csv")
RNA_Lip_DESeq_results_autosomes_F <- RNA_Lip_sex_stratified[RNA_Lip_sex_stratified$comparison == "F_PEvsF_Cont", ]
RNA_Lip_DESeq_results_autosomes_M <- RNA_Lip_sex_stratified[RNA_Lip_sex_stratified$comparison == "M_PEvsM_Cont", ]
RNA_WG_sex_stratified <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_WG_rerun/RNA_WG_DESeq_results_autosomes_IT_all.csv")
RNA_WG_DESeq_results_autosomes_F <- RNA_WG_sex_stratified[RNA_WG_sex_stratified$comparison == "F_PEvsF_Cont", ]
RNA_WG_DESeq_results_autosomes_M <- RNA_WG_sex_stratified[RNA_WG_sex_stratified$comparison == "M_PEvsM_Cont", ]



Swish_isoform_combined <- Swish_isoform_combined[Swish_isoform_combined$Q_Value <= 0.05, c("Transcript_ID", "Gene_ID.value", "Log2_FoldChg", "Q_Value")]

#Make Differential Isoform (Swish) and Isoform Switch (IsoformSwitchAnalyzeR) merged table WG
Swish_IsoformSwitch_WG <- Isoform_Switch_WG
Swish_IsoformSwitch_WG$isoformUpregulated_Sig_Swish_Log2FC <- Swish_isoform_combined$Log2_FoldChg[match(Swish_IsoformSwitch_WG$isoformUpregulated, Swish_isoform_combined$Transcript_ID)]
Swish_IsoformSwitch_WG$isoformDownregulated_Sig_Swish_Log2FC <- Swish_isoform_combined$Log2_FoldChg[match(Swish_IsoformSwitch_WG$isoformDownregulated, Swish_isoform_combined$Transcript_ID)]
Swish_IsoformSwitch_WG$Gene_ID.value <- Swish_isoform_combined$Gene_ID.value[match(Swish_IsoformSwitch_WG$isoformDownregulated, Swish_isoform_combined$Transcript_ID)]
write.csv(Swish_IsoformSwitch_WG, file = "Swish_IsoformSwitch_WG.csv")

#Make Differential Isoform (Swish) and Isoform Switch (IsoformSwitchAnalyzeR) + DGE merged table WG
DTE_Switch_DGE_WG <- merge(Swish_IsoformSwitch_WG, RNA_WG_DESeq_results_autosomes_combined_sex, by.x = "Gene_ID.value", by.y = "Row.names", all.x = TRUE)
write.csv(DTE_Switch_DGE_WG,"DTE_Switch_DGE_WG.csv")

#Make Differential Isoform (Swish) and Isoform Switch (IsoformSwitchAnalyzeR) + DGE + rMATS merged table WG
DTE_Switch_DGE_rMATS_WG <- merge(DTE_Switch_DGE_WG, rMATS_results_WG, by.x = "Gene_ID.value", by.y = "GeneID", all.x = TRUE)
write.csv(DTE_Switch_DGE_rMATS_WG,"rMATS_DTE_Switch_DGE_WG.csv") #0



#Make formation of tables into a function

#Make single table with all data

Isoform_Table <- function (SwitchList, file_prefix, Swish_name, DESeq2_name, rMATS_name){
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



Swish_name <- Swish_name[Swish_name$Q_Value <= 0.05, c("Transcript_ID", "Gene_ID.value", "Log2_FoldChg", "Q_Value")]

#Make Differential Isoform (Swish) and Isoform Switch (IsoformSwitchAnalyzeR) merged table WG
Swish_top_switches <- top_switches
Swish_top_switches$isoformUpregulated_Sig_Swish_Log2FC <- Swish_name$Log2_FoldChg[match(Swish_top_switches$isoformUpregulated, Swish_name$Transcript_ID)]
Swish_top_switches$isoformDownregulated_Sig_Swish_Log2FC <- Swish_name$Log2_FoldChg[match(Swish_top_switches$isoformDownregulated, Swish_name$Transcript_ID)]
# Swish_top_switches$Gene_ID.value <- Swish_name$Gene_ID.value[match(Swish_top_switches$isoformDownregulated, Swish_name$Transcript_ID)]
Swish_top_switches$Gene_ID.value <- gtf_data$gene_id[match(sub(";.*", "", Swish_top_switches$isoformUpregulated), gtf_data$transcript_id)]

filename_DTE_Switch<- paste0(file_prefix, "_Swish_IsoformSwitch.csv")

# write.csv(Swish_top_switches, filename_DTE_Switch)

#Make Differential Isoform (Swish) and Isoform Switch (IsoformSwitchAnalyzeR) + DGE + rMATS merged table WG
DTE_Switch_DGE <- merge(Swish_top_switches, DESeq2_name, by.x = "Gene_ID.value", by.y = "Row.names", all.x = TRUE)
DTE_Switch_DGE$GeneID <- sub("\\..*", "", DTE_Switch_DGE$Gene_ID.value)
DTE_Switch_DGE_rMATS <- merge(DTE_Switch_DGE, rMATS_name, by.x = "Gene_ID.value", by.y = "GeneID", all.x = TRUE)

filename_DTE_Switch_DGE_rMATS <- paste0(file_prefix, "_Swish_IsoformSwitch.csv")
write.csv(DTE_Switch_DGE_rMATS, filename_DTE_Switch_DGE_rMATS) #0


#Make Comparisons only between Swish, rMATS and DESeq2
Swish_rMATS <-  merge(Swish_name, DESeq2_name, by.x = "Gene_ID.value", by.y = "Row.names", all.x = TRUE)
Swish_rMATS$GeneID <- sub("\\..*", "", Swish_rMATS$Gene_ID.value)
Swish_rMATS_DEG <- merge(Swish_rMATS, rMATS_name, by.x = "Gene_ID.value", by.y = "GeneID", all.x = TRUE)

filename_Swish_rMATS_DEG <- paste0(file_prefix, "_Swish_rMATS_DEG.csv")
write.csv(Swish_rMATS_DEG, filename_Swish_rMATS_DEG) 


#Overview Volcano Plot

## Each dot is an Isoform that has been switched (and sig and non-sig is indicated by low dIF)

vol_plot_filename <- paste0("./", file_prefix, "Isoform_Switch_Vol.png")

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

switch_v_gene_filename <- paste0("./", file_prefix, "Isoform_Switch_vs_Gene_Vol.png")
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

dir.create("./WG_combined_sex_Isoform_Tables")
setwd("./WG_combined_sex_Isoform_Tables")
Isoform_Table(exampleSwitchListAnalyzed, "WG_combined_sex", Swish_isoform_combined, RNA_WG_DESeq_results_autosomes_combined_sex, rMATS_results_WG)



dir.create("./Lip_combined_sex_Isoform_Tables")
setwd("./Lip_combined_sex_Isoform_Tables")

Isoform_Table(exampleSwitchListAnalyzed_Lip, "Lip_combined_sex", Swish_Lip_isoform_combined, RNA_Lip_DESeq_results_autosomes_combined_sex, rMATS_results_WG)













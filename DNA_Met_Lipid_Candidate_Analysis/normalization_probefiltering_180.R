BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19", update = F)
BiocManager::install("IlluminaHumanMethylation450kmanifest")
BiocManager::install("minfi")
BiocManager::install("magrittr")
BiocManager::install("wateRmelon")
library(minfi)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(magrittr)
library(wateRmelon)

cd /workspace/lab/wilsonslab/datalake-wilsonslab/PE_Lipid_Meta-analysis/2025_DNAm_Lipid_Candidate/Candidate_Diff_DNAm_Analysis_2025/normalization_probefiltering/

#Load Base Directories (each directory should have all red and green IDATs + one sample sheet as a .csv)
baseDir98224 <- ("/workspace/lab/wilsonslab/datalake-wilsonslab/PE_Lipid_Meta-analysis/2025_DNAm_Lipid_Candidate/raw_data/GSE98224")
baseDir100197 <- ("/workspace/lab/wilsonslab/datalake-wilsonslab/PE_Lipid_Meta-analysis/2025_DNAm_Lipid_Candidate/raw_data/GSE100197")
baseDir125605 <- ("/workspace/lab/wilsonslab/datalake-wilsonslab/PE_Lipid_Meta-analysis/2025_DNAm_Lipid_Candidate/raw_data/GSE125605")
baseDir75196 <- ("/workspace/lab/wilsonslab/datalake-wilsonslab/PE_Lipid_Meta-analysis/2025_DNAm_Lipid_Candidate/raw_data/GSE75196")

#Create targets 
targets98224 <- read.metharray.sheet(baseDir98224)
targets100197 <- read.metharray.sheet(baseDir100197)
targets125605 <- read.metharray.sheet(baseDir125605)
targets75196 <- read.metharray.sheet(baseDir75196)

#Exclude samples 
edit_targets98224 <- targets98224
edit_targets100197 <- subset(targets100197, Excluded !="Excluded")
edit_targets125605 <- subset(targets125605, Excluded !="Excluded")
edit_targets75196 <- targets75196
# Only excluded samples in GSE_100197 where some placentas were IUGR (n=11), pre-term (n=24), or technical replicates (n=8) and GSE_125605 where there was no gestational age (n=1) - based on metadata information

#Build RG sets
 
RGset100197 <- read.metharray.exp(targets = edit_targets100197, verbose = TRUE)
RGset125605 <- read.metharray.exp(targets = edit_targets125605, verbose = TRUE)
RGset75196 <- read.metharray.exp(targets = edit_targets75196, verbose = TRUE)
RGset98224 <- read.metharray.exp(targets = edit_targets98224, verbose = TRUE)

#Combine the RG sets (can only do this two at a time)
combo_100197_125605 <- combineArrays (RGset100197,RGset125605, outType = c("IlluminaHumanMethylation450k"), verbose = TRUE)
combo_100197_125605_75196 <- combineArrays (combo_100197_125605,RGset75196, outType = c("IlluminaHumanMethylation450k"), verbose = TRUE)
combined_RGset <-  combineArrays (combo_100197_125605_75196,RGset98224, outType = c("IlluminaHumanMethylation450k"), verbose = TRUE)

#Check number of probes (type I and type II probes added together = 485,512 for 450K)
combined_manifest <- getManifest(combined_RGset)
print(combined_manifest)
    IlluminaMethylationManifest object
    Annotation
      array: IlluminaHumanMethylation450k
    Number of type I probes: 135476
    Number of type II probes: 350036
    Number of control probes: 850
    Number of SNP type I probes: 25
    Number of SNP type II probes: 40

#Normalization (adjFunnorm), Background and dye bias correction with noob, mapping to genome, quantile extraction, normalization
saveRDS(combined_RGset, file = "placmet_noNorm.RDS")
placmet_adjFunnorm <- adjustedFunnorm(combined_RGset, verbose = TRUE)
saveRDS(placmet_adjFunnorm, file = "placmet_adjFunnorm.RDS")

#Adapted filtering step to keep both autosomal and XY probes

#Filtering bad probes (bad detection p vals or missing betas)
detp <- detectionP(combined_RGset)
number_bad_P_before_Fstrat <- print(sum(rowSums(detp)>=(ncol(combined_RGset))*0.05)) #422
write.csv(detp, "detp_table.csv")

## Annotate X and Y probes
annotation_price38 <- read.csv("/workspace/lab/wilsonslab/datalake-wilsonslab/PE_Lipid_Meta-analysis/2025_DNAm_Lipid_Candidate/Candidate_DNAm_Annotations/full_annotation_region_vs_price38.csv")
chrXprobes <- subset(annotation_price38, annotation_price38$chr == "chrX") #8294
chrYprobes <- subset(annotation_price38, annotation_price38$chr == "chrY") #185

# probeInfo <- as.data.frame(cbind(IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Locations, IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Other, IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Manifest)) 
# probeInfo$probeID <- rownames(probeInfo)
# chrXprobes_pI <- subset(probeInfo, probeInfo$chr == "chrX") #19090
# chrYprobes_pI <- subset(probeInfo, probeInfo$chr == "chrY") #537

#metadata <- read.csv("/workspace/lab/wilsonslab/datalake-wilsonslab/PE_Lipid_Meta-analysis/2025_DNAm_Lipid_Candidate/Candidate_Diff_DNAm_Analysis_2025/Metadata_Sheet_lipid_preeclampsia_excluded_removed.csv")
metadata <- read.csv("/workspace/lab/wilsonslab/datalake-wilsonslab/PE_Lipid_Meta-analysis/2025_DNAm_Lipid_Candidate/Candidate_Diff_DNAm_Analysis_2025/totalsamples_180/metadata.csv")


males <- subset(metadata, metadata$Fetal_Sex == "M")
females <- subset(metadata, metadata$Fetal_Sex == "F")
## For females, set detp in Y chromosomes to 0 (these probes do not bind to anything in female samples)
detp[rownames(detp) %in% chrYprobes$CpG_ID, females$Sample_Name] <- 0
head(detp[rownames(detp) %in% chrYprobes$CpG_ID, females$Sample_Name]) #Make sure all 0
                        GSM1944959_9376561070_R05C01 GSM1944960_9376561070_R06C01
         cg00212031                            0                            0
         cg00213748                            0                            0
         cg00214611                            0                            0
         cg00455876                            0                            0
         cg01707559                            0                            0
         cg02004872                            0                            0
                    GSM1944966_9376561070_R06C02
         cg00212031                            0                                     
         cg00213748                            0
         cg00214611                            0
         cg00455876                            0
         cg01707559                            0
         cg02004872                            0

# bad probes have detection p-value > 0.01
bad_detp <- detp > 0.01
#number of bad probes in >= 5% of samples
number_bad_detp <- print(sum(rowSums(bad_detp)>=(ncol(combined_RGset))*0.05)) #745
# missing betas >= 5% of samples 
avgbeta <- getBeta(combined_RGset)
bad_beta <- is.na(avgbeta)
number_bad_beta <- print(sum(rowSums(bad_beta)>=(ncol(combined_RGset))*0.05)) #17

# remove probes with bad p-values or missing betas >= 5% samples
badProbes <- rowSums(bad_detp)>=(ncol(combined_RGset))*0.05 | rowSums(bad_beta)>=(ncol(combined_RGset))*0.05 
badProbes1 <- as.data.frame(badProbes) 
#rowSums function below didn't work on badProbes because it's a vector so turned it into a dataframe
number_badprobes <- print(sum(rowSums(badProbes1, na.rm = TRUE))) #746
placmet_adjFunnorm_BPfilt <- placmet_adjFunnorm[!badProbes,]

# Remove SNP probes
# download full table from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL16304
price_anno <- read.delim("/workspace/lab/wilsonslab/eyerk/2025_Lipid_Candidate_GSE_Info/GSE_annotations/GPL16304-47833_no_legend.tsv",header=TRUE)
price_anno_SNP <- subset(price_anno, price_anno$n_target.CpG.SNP>0)
placmet_adjFunnorm_SNPremoved <- placmet_adjFunnorm_BPfilt[!rownames(placmet_adjFunnorm_BPfilt) %in% price_anno_SNP$ID,]

# Remove Cross-Hybridizing probes
price_anno_Hyb <- subset(price_anno, price_anno$XY_Hits == "XY_YES" | price_anno$Autosomal_Hits == "A_YES")
placmet_adjFunnorm_HybRemoved <- placmet_adjFunnorm_SNPremoved[!rownames(placmet_adjFunnorm_SNPremoved) %in% price_anno_Hyb$ID,]

# Remove Non-variable placental probes, #downloaded from https://github.com/redgar598/Tissue_Nonvariable_450K_CpGs
library(RCurl)
x <- getURL("https://raw.githubusercontent.com/redgar598/Tissue_Invariable_450K_CpGs/master/Invariant_Placenta_CpGs.csv")
write.csv (x,"Invariant_Placenta_CpGs.csv")
nonvarplac_anno <- read.csv("./Invariant_Placenta_CpGs.csv", sep = ",") 
placmet_adjFunnorm_nonvarremoved <- placmet_adjFunnorm_HybRemoved[!rownames(placmet_adjFunnorm_HybRemoved) %in% nonvarplac_anno$CpG,]
placmet_adjFunnorm_allfiltered <- placmet_adjFunnorm_nonvarremoved

saveRDS(placmet_adjFunnorm_allfiltered, "placmet_adjFunnorm_allfiltered_180.rds")


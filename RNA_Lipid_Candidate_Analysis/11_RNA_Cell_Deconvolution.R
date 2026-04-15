#Cellular Deconvolution via MuSiC2
#Goal: Estimate cell type proportions from bulk RNA-seq using a multi-subject single-cell reference to account for inter-individual variability.

install.packages(c("cli", "dplyr", "forcats", "ggplot2", "lifecycle", "patchwork", "purrr", "rlang", "scales", "stats", "stringr", "utils", "tidyr"))
install.packages("ggstats")
install.packages("GGally")
install.packages('locfdr')
install.packages(c('nnls', 'corpcor'))
install.packages("EpiDISH_2.26.0.tar.gz")
install.packages("TOAST_1.25.0.tar.gz")
BiocManager::install(c("rtracklayer", "ggpubr", "tidyverse", "planet", "minfi", "EpiDISH"))

devtools::install_github('xuranw/MuSiC')
library(MuSiC)
library(rtracklayer)
library(data.table)
library(dplyr)
library(tidyr)
library(Biobase)
library(SingleCellExperiment)
library(lme4) #Version 1.1.35.4
library(lmerTest) #Version 3.1.3

#Download reference placental scRNA-seq dataset in Bash 
wget https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE182381&format=file&file=GSE182381_reference_sample.txt.gz
gunzip GSE182381_reference_sample.txt.gz

#Load RNA Read Counts
read_counts <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/read_counts/read_count_nohead.csv")
RNA_metadata <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/metadata/RNA_Lipid_Candidate_Metadata.csv")
reference_scRNAseq <- read.delim("./GSE182381_reference_sample.txt/GSE182381_reference_sample.txt")

#alter Read Counts file to be in correct format (1st column is titled 'gene names' and row values are gene name values, subsequent column names are sample names)
gtf_file <- "/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/genome_mapping/gencode.v48.primary_assembly.annotation.gtf"
gtf_data <- import(gtf_file)
gtf_df <- as.data.frame(gtf_data)
unique_gtf_df <- gtf_df[!duplicated(gtf_df$gene_id),]

metadata_keep <- RNA_metadata[RNA_metadata$exclude == "keep",]

colnames(read_counts) <- read_counts [1,]
read_counts <- read_counts [-(1), ]
colnames(read_counts) <- sub("../genome_mapping/markeddup_BAMs/", "", colnames(read_counts))
colnames(read_counts) <- sub("_markdup.bam", "", colnames(read_counts))

read_counts_samples <- read_counts[,metadata_keep$Run %in% colnames(read_counts[7:243]) ]
read_counts_genename <- left_join(read_counts_samples, unique_gtf_df[, c("gene_id", "gene_name")], by = c("Geneid" = "gene_id"))
read_counts_final <- read_counts_genename [,-(1:6)]

# Since reference file counts only expressed in gene name, need to merge counts between ENSEMBLIDs with same gene name
sum(duplicated(read_counts_final$gene_name))
numeric <- colnames(read_counts_final)
numeric <- numeric[1:122]
read_counts_final[numeric] <- lapply(read_counts_final[numeric], as.numeric)

read_counts_final <- read_counts_final %>%
  group_by(across(all_of("gene_name"))) %>%
  summarise(across(all_of(numeric), sum), 
.groups = "drop")
read_counts_final_df <- as.data.frame(read_counts_final)

row.names(read_counts_final_df) <- read_counts_final_df$gene_name
read_counts_final_df <- read_counts_final_df[,(2:123)]
write.csv(read_counts_final_df, "read_counts_final.csv")


# MuSiC Input Data Setup

## Set-up Reference Dataset for Analysis
gene_exprs_mtx <- reference_scRNAseq
row.names(gene_exprs_mtx) <- gene_exprs_mtx$GeneSymbol 
gene_exprs_mtx <- gene_exprs_mtx [,-1]
gene_exprs_mtx_all <- as.matrix(gene_exprs_mtx)
gene_exprs_mtx_condensed <- gene_exprs_mtx_all[rownames(gene_exprs_mtx_all) %in% rownames(read_counts_final_df), ]

pheno_matrix <- data.frame(row.names = colnames(gene_exprs_mtx_condensed))
pheno_matrix$cell_id <- row.names(pheno_matrix)
pheno_matrix$cell_id <- sub("\\.[0-9]+$", "", pheno_matrix$cell_id)
pheno_matrix$subject_name <- "control" #Disease Group
pheno_matrix$sample_id <- "PseudoSample_1" #need to create multiple fake sample_ids for algorithm to work, randomly group into 3 psuedogroups
pheno_matrix$sample_id[13499:26996] <- "PseudoSample_2" 
pheno_matrix$sample_id[26997:40494] <- "PseudoSample_3" 

sc_sce <- SingleCellExperiment(list(counts = gene_exprs_mtx_condensed), colData = pheno_matrix)


## Set-up Bulk Dataset for Analysis
metadata_keep_eset <- metadata_keep
row.names(metadata_keep_eset) <- metadata_keep_eset$Run
metadata_keep_eset_control <- metadata_keep_eset[metadata_keep_eset$disease_group == "Control", ]
metadata_keep_eset_PE <- metadata_keep_eset[metadata_keep_eset$disease_group == "PE", ]

read_counts_final_control <- read_counts_final_df[rownames(read_counts_final_df) %in% rownames(gene_exprs_mtx_condensed) , colnames(read_counts_final_df) %in% metadata_keep_eset_control$Run]
read_counts_final_PE <- read_counts_final_df[rownames(read_counts_final_df) %in% rownames(gene_exprs_mtx_condensed), colnames(read_counts_final_df) %in% metadata_keep_eset_PE$Run]


bulk_metadata <- data.frame(
    labelDescription = c("X", "geo_accession", "title", "source_name_ch1", "molecule_ch1", "platform_id", "instrument_model", "library_strategy", "GSE_number", "age.of_mother_.years..ch1", "term.ch1", "sga.ch1", "rop.ch1", "Run", "BioProject", "BioSample", "SRA.Study", "LibraryLayout", "disease_group", "fetal_sex", "gestational_age_weeks_days", "maternal_ethnicity", "exclude", "predicted_fetal_sex"),
    row.names = c("X", "geo_accession", "title", "source_name_ch1", "molecule_ch1", "platform_id", "instrument_model", "library_strategy", "GSE_number", "age.of_mother_.years..ch1", "term.ch1", "sga.ch1", "rop.ch1", "Run", "BioProject", "BioSample", "SRA.Study", "LibraryLayout", "disease_group", "fetal_sex", "gestational_age_weeks_days", "maternal_ethnicity", "exclude", "predicted_fetal_sex")
)

plac_bulk_eset_control <- ExpressionSet(assayData = data.matrix(read_counts_final_control), phenoData = new("AnnotatedDataFrame", data = metadata_keep_eset_control, varMetadata = bulk_metadata))
plac_bulk_eset_PE <- ExpressionSet(assayData = data.matrix(read_counts_final_PE), phenoData = new("AnnotatedDataFrame", data = metadata_keep_eset_PE, varMetadata = bulk_metadata))

bulk_control_mtx <- exprs(plac_bulk_eset_control)
bulk_PE_mtx <- exprs(plac_bulk_eset_PE)

# Running Cell Deconvolution
Est_prop_plac <-  music2_prop_t_statistics(bulk.control.mtx = bulk_control_mtx, bulk.case.mtx = bulk_PE_mtx, sc.sce = sc_sce, clusters = 'cell_id',
                               samples = 'sample_id', select.ct = c("Fetal.Mesenchymal.Stem.Cells", "Fetal.CD14..Monocytes", "Fetal.CD8..Activated.T.Cells", "Fetal.Naive.CD4..T.Cells", "Fetal.Naive.CD8..T.Cells", "Fetal.Natural.Killer.T.Cells", "Fetal.B.Cells", "Fetal.GZMK..Natural.Killer", "Fetal.Memory.CD4..T.Cells", "Fetal.Hofbauer.Cells", "Fetal.Plasmacytoid.Dendritic.Cells", "Fetal.GZMB..Natural.Killer", "Fetal.Endothelial.Cells", "Fetal.Syncytiotrophoblast", "Fetal.Fibroblasts", "Fetal.Cytotrophoblasts", "Fetal.Proliferative.Cytotrophoblasts", "Fetal.Nucleated.Red.Blood.Cells", "Maternal.CD8..Activated.T.Cells", "Maternal.Naive.CD4..T.Cells", "Maternal.FCGR3A..Monocytes", "Maternal.CD14..Monocytes", "Maternal.Natural.Killer.Cells", "Maternal.B.Cells", "Maternal.Plasma.Cells", "Maternal.Naive.CD8..T.Cells", "Fetal.Extravillous.Trophoblasts"
))



Final_est_plac_prop <- Est_prop_plac$Est.prop
Final_est_plac_prop <- as.data.frame(Final_est_plac_prop)
write.csv(Final_est_plac_prop, file = "Final_est_plac_prop.csv")

m_CONT_F <- metadata_keep_eset_control[metadata_keep_eset_control$predicted_fetal_sex == "F",]
m_CONT_M <- metadata_keep_eset_control[metadata_keep_eset_control$predicted_fetal_sex == "M",]
m_PE_F <- metadata_keep_eset_PE[metadata_keep_eset_PE$predicted_fetal_sex == "F",]
m_PE_M <- metadata_keep_eset_PE[metadata_keep_eset_PE$predicted_fetal_sex == "M",]

estF_CONT_F <- Final_est_plac_prop[rownames(Final_est_plac_prop) %in% rownames(m_CONT_F),]
estF_CONT_M <- Final_est_plac_prop[rownames(Final_est_plac_prop) %in% rownames(m_CONT_M),]
estF_PE_F <- Final_est_plac_prop[rownames(Final_est_plac_prop) %in% rownames(m_PE_F),]
estF_PE_M <- Final_est_plac_prop[rownames(Final_est_plac_prop) %in% rownames(m_PE_M),]

estF_PE_F$group <- 'PE_F'
estF_CONT_F$group <- 'CONT_F'
estF_PE_M$group <- 'PE_M'
estF_CONT_M$group <- 'CONT_M'

estF <- rbind(estF_PE_F, estF_CONT_F, estF_PE_M, estF_CONT_M)
# estF <- as.numeric(estF[1:27])



estF_diseaseandsex <- estF %>%
  mutate(disease_group = case_when( #Creating the disease_group variable column
    group == "PE_F" | group == "PE_M" ~ "PE",
    group == "CONT_F" | group == "CONT_M" ~ "CONT",
    .default = NA)) %>%
  mutate(sex = case_when( #Creating the sex variable column
    group == "PE_F" | group == "CONT_F" ~ "F",
    group == "PE_M" | group == "CONT_M" ~ "M",
    .default = NA))
estF_diseaseandsex$disease_group <- as.factor(estF_diseaseandsex$disease_group)
estF_diseaseandsex$sex <- as.factor(estF_diseaseandsex$sex)
estF_diseaseandsex$Sample_Name <- rownames(estF_diseaseandsex)
metadata_keep_eset$Sample_Name <- rownames(metadata_keep_eset)
estF_dat <- merge(estF_diseaseandsex, metadata_keep_eset, by = "Sample_Name")
estF_dat$predicted_fetal_sex <- as.factor(estF_dat$predicted_fetal_sex)
names(estF_dat)[names(estF_dat) == "disease_group.x"] <- "disease_group"

#Filter-out Negligible Cell Proportions
all_cell_types <- c("Fetal.Mesenchymal.Stem.Cells", "Fetal.CD14..Monocytes", "Fetal.CD8..Activated.T.Cells", "Fetal.Naive.CD4..T.Cells", "Fetal.Naive.CD8..T.Cells", "Fetal.Natural.Killer.T.Cells", "Fetal.B.Cells", "Fetal.GZMK..Natural.Killer", "Fetal.Memory.CD4..T.Cells", "Fetal.Hofbauer.Cells", "Fetal.Plasmacytoid.Dendritic.Cells", "Fetal.GZMB..Natural.Killer", "Fetal.Endothelial.Cells", "Fetal.Syncytiotrophoblast", "Fetal.Fibroblasts", "Fetal.Cytotrophoblasts", "Fetal.Proliferative.Cytotrophoblasts", "Fetal.Nucleated.Red.Blood.Cells", "Maternal.CD8..Activated.T.Cells", "Maternal.Naive.CD4..T.Cells", "Maternal.FCGR3A..Monocytes", "Maternal.CD14..Monocytes", "Maternal.Natural.Killer.Cells", "Maternal.B.Cells", "Maternal.Plasma.Cells", "Maternal.Naive.CD8..T.Cells", "Fetal.Extravillous.Trophoblasts")
cell_means <- colMeans(estF_dat[, all_cell_types]) 
cell_types <- names(cell_means[cell_means >= 0.0001])
colMeans(estF_dat[, cell_types]) 

#Linear Model Comparing Cell Proportions Between Groups 
#Males: Comparing PE sample cell proportions to Control sample cell proportions
male_cell_results <- data.frame(
  Cell_Type = character(),
  Beta_PathGroup = numeric(),
  P_PathGroup = numeric(),
  stringsAsFactors = FALSE
)
male_data <- subset(estF_dat, predicted_fetal_sex == "M")
for(cell_type in cell_types) {
  cat("Analyzing cell types:", cell_type, "\n")
  
  #Fit linear model (comparing conception type)
  formula <- as.formula(paste(cell_type, "~ disease_group + (1 | GSE_number)"))
  lmm <- lmer(formula, data = male_data)
  
  #model summary 
  model_summary <- summary(lmm)
  
  #Extract fixed effects 
  fixed_effects <-  model_summary$coefficients 
  
  #extract results for conception type 
  beta_disease_group <- fixed_effects["disease_groupPE", "Estimate"]
  p_disease_group <- fixed_effects["disease_groupPE", "Pr(>|t|)"]
  
  # Append results to the data frame
  male_cell_results <- rbind(male_cell_results, data.frame(
    Cell_Type = cell_type,
    Beta_PathGroup = beta_disease_group,
    P_PathGroup = p_disease_group
  ))
}
write.csv(male_cell_results, "male_stratified_celltype_results.csv", row.names = FALSE)
#Females: Comparing PE sample cell proportions to Control sample cell proportions
female_cell_results <- data.frame(
  Cell_Type = character(),
  Beta_PathGroup = numeric(),
  P_PathGroup = numeric(),
  stringsAsFactors = FALSE
)
female_data <- subset(estF_dat, predicted_fetal_sex == "F")
for(cell_type in cell_types) {
  cat("Analyzing cell types:", cell_type, "\n")
  
  #Fit linear model (comparing conception type)
  formula <- as.formula(paste(cell_type, "~ disease_group + (1 | GSE_number)"))
  lmm <- lmer(formula, data = female_data)
  
  #model summary 
  model_summary <- summary(lmm)
  
  #Extract fixed effects 
  fixed_effects <-  model_summary$coefficients 
  
  #extract results for conception type 
  beta_disease_group <- fixed_effects["disease_groupPE", "Estimate"]
  p_disease_group <- fixed_effects["disease_groupPE", "Pr(>|t|)"]
  
  # Append results to the data frame
  female_cell_results <- rbind(female_cell_results, data.frame(
    Cell_Type = cell_type,
    Beta_PathGroup = beta_disease_group,
    P_PathGroup = p_disease_group
  ))
}
write.csv(female_cell_results, "female_stratified_celltype_results.csv", row.names = FALSE)
#PE only Samples: Comparing Male and Female Samples
PE_cell_results <- data.frame(
  Cell_Type = character(),
  Beta_Sex = numeric(),
  P_Sex = numeric(),
  stringsAsFactors = FALSE
)
PE_data <- subset(estF_dat, disease_group == "PE")
for(cell_type in cell_types) {
  cat("Analyzing cell types:", cell_type, "\n")
  
  #Fit linear model (comparing conception type)
  formula <- as.formula(paste(cell_type, "~ predicted_fetal_sex + (1 | GSE_number)"))
  lmm <- lmer(formula, data = PE_data)
  
  #model summary 
  model_summary <- summary(lmm)
  
  #Extract fixed effects 
  fixed_effects <-  model_summary$coefficients 
  
  #extract results for conception type 
  beta_sex <- fixed_effects["predicted_fetal_sexM", "Estimate"]
  p_sex <- fixed_effects["predicted_fetal_sexM", "Pr(>|t|)"]
  
  # Append results to the data frame
  PE_cell_results <- rbind(PE_cell_results, data.frame(
    Cell_Type = cell_type,
    Beta_Sex = beta_sex,
    P_Sex = p_sex
  ))
}
write.csv(PE_cell_results, "PE_stratified_celltype_results.csv", row.names = FALSE)
#Control only Samples: Comparing Male and Female Samples
CONT_cell_results <- data.frame(
  Cell_Type = character(),
  Beta_Sex = numeric(),
  P_Sex = numeric(),
  stringsAsFactors = FALSE
)
CONT_data <- subset(estF_dat, disease_group == "CONT")
for(cell_type in cell_types) {
  cat("Analyzing cell types:", cell_type, "\n")
  
  #Fit linear model (comparing conception type)
  formula <- as.formula(paste(cell_type, "~ predicted_fetal_sex + (1 | GSE_number)"))
  lmm <- lmer(formula, data = CONT_data)
  
  #model summary 
  model_summary <- summary(lmm)
  
  #Extract fixed effects 
  fixed_effects <-  model_summary$coefficients 
  
  #extract results for conception type 
  beta_sex <- fixed_effects["predicted_fetal_sexM", "Estimate"]
  p_sex <- fixed_effects["predicted_fetal_sexM", "Pr(>|t|)"]
  
  # Append results to the data frame
  CONT_cell_results <- rbind(CONT_cell_results, data.frame(
    Cell_Type = cell_type,
    Beta_Sex = beta_sex,
    P_Sex = p_sex
  ))
}
write.csv(CONT_cell_results, "CONT_stratified_celltype_results.csv", row.names = FALSE)




#Mixed Effects Model - No stratification of sex or conception type done prior
# Initialize a results data frame
cell_unstratified_results <- data.frame(
  Cell_Type = character(),
  Beta_Sex = numeric(),
  P_Sex = numeric(),
  Beta_PathGroup = numeric(),
  P_PathGroup = numeric(),
  stringsAsFactors = FALSE
)
# Loop through each cell type
for (cell_type in cell_types) {
  cat("Analyzing cell type:", cell_type, "\n")
  
  # Fit the linear mixed-effects model
  formula <- as.formula(paste(cell_type, "~ predicted_fetal_sex + disease_group  + (1 | GSE_number)"))
  lmm <- lmer(formula, data = estF_dat)
  
  # Summarize the model
  model_summary <- summary(lmm)
  
  # Extract fixed effects
  fixed_effects <- model_summary$coefficients
  
  # Extract results for sex and conception type
  beta_sex <- fixed_effects["predicted_fetal_sexM", "Estimate"]
  p_sex <- fixed_effects["predicted_fetal_sexM", "Pr(>|t|)"]
  
  beta_disease_group <- fixed_effects["disease_groupPE", "Estimate"]
  p_disease_group <- fixed_effects["disease_groupPE", "Pr(>|t|)"]
  
  # Append results to the data frame
  cell_unstratified_results <- rbind(cell_unstratified_results, data.frame(
    Cell_Type = cell_type,
    Beta_Sex = beta_sex,
    P_Sex = p_sex,
    Beta_PathGroup = beta_disease_group,
    P_PathGroup = p_disease_group
  ))
}
write.csv(cell_unstratified_results, "LMM_CellType_Results.csv", row.names = FALSE)

#Plotting 
#Reshape data into long format for ggplot2
#Plot separated by study 
long_data <- estF_dat %>%
  pivot_longer(
    cols = c(Fetal.Syncytiotrophoblast, Fetal.Cytotrophoblasts, Fetal.Proliferative.Cytotrophoblasts, Fetal.Extravillous.Trophoblasts, Fetal.Fibroblasts, Fetal.Mesenchymal.Stem.Cells, Fetal.Hofbauer.Cells, Fetal.Endothelial.Cells, Fetal.Nucleated.Red.Blood.Cells, Fetal.GZMK..Natural.Killer, Fetal.GZMB..Natural.Killer),
    names_to = "Cell_Type",
    values_to = "Proportion"
  )
png("celltype.png", width=4000, height=1500, res=200)
ggplot(long_data, aes(x = Cell_Type, y = Proportion, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(group = group), width = 0.2, alpha = 0.1, size = 1) +
  facet_wrap(~ GSE_number, scales = "free_y") +
  labs(
    title = "Boxplot of Cell Type Proportions by Group and Sex - RNA-Seq",
    x = "Group",
    y = "Proportion",
    fill = "Sex"
  ) +
  theme_minimal() +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
dev.off()
#plot with all samples in one figure
png("celltype_oneplot.png",  height = 7.5, width = 10, units = "in", res = 750)
ggplot(long_data, aes(x = Cell_Type, y = Proportion, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +
  geom_jitter(aes(group = group), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75),alpha = 0.1, size = 1) +
  labs(
    title = "Boxplot of Cell Type Proportions by Group - RNA-Seq",
    x = "Cell Type",
    y = "Proportion",
    fill = "Group"
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )
dev.off()








































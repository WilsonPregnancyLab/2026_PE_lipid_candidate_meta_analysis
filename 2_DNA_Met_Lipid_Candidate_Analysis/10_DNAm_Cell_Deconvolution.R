BiocManager::install(c("lme4", "lmerTest", "ggsignf", "ggpubr"))

library(ggplot2) #Version 3.5.1
library(dplyr) #Version 1.1.4
library(tidyr) #Version 1.3.1
library(lme4) #Version 1.1.35.4
library(lmerTest) #Version 3.1.3
library(tidyr) #Version 1.3.1


library(ggsignif) #version 0.6.4
library(ggpubr) #version 1.0.0
library(tidyverse) #version 2.0.0
library(planet) #version 1.20.0
library(minfi) #version 1.58.0
library(EpiDISH) #version 2.28.0
data("plBetas")
data("plCellCpGsThird")

#Creating groups for sex stratification
metadata <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_DNAm_Lipid_Candidate/Candidate_Diff_DNAm_Analysis_2025/Metadata_Sheet_lipid_preeclampsia_excluded_removed.csv")
metadata$pathology_group <- as.factor(metadata$pathology_group)
metadata$Fetal_Sex <- as.factor(metadata$Fetal_Sex)
metadata$Sentrix_ID <- as.factor(metadata$Sentrix_ID)
metadata$Sentrix_Position <- as.factor(metadata$Sentrix_Position)
metadata$GSE_number <- as.factor(metadata$GSE_number)
metadata$geo_accession <- as.factor(metadata$geo_accession)
metadata$gestational_age <- as.factor(metadata$gestational_age)

males <- subset(metadata, Fetal_Sex == "M")
males$pathology_group <- as.factor(males$pathology_group)
males$Fetal_Sex <- as.factor(males$Fetal_Sex)
males$Sentrix_ID <- as.factor(males$Sentrix_ID)
males$Sentrix_ID <- droplevels(males$Sentrix_ID)
males$Sentrix_Position <- as.factor(males$Sentrix_Position)
males$GSE_number <- as.factor(males$GSE_number)

females <- subset(metadata, Fetal_Sex == "F")
females$pathology_group <- as.factor(females$pathology_group)
females$Fetal_Sex <- as.factor(females$Fetal_Sex)
females$Sentrix_ID <- as.factor(females$Sentrix_ID)
females$Sentrix_ID <- droplevels(females$Sentrix_ID)
females$Sentrix_Position <- as.factor(females$Sentrix_Position)
females$GSE_number <- as.factor(females$GSE_number)

# Read the filtered and normalized RG sets
placmet_adjFunnorm_allfiltered <- readRDS("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_DNAm_Lipid_Candidate/Candidate_Diff_DNAm_Analysis_2025/placmet_adjFunnorm_allfiltered.rds") #dim 424750 172
placmet_adjFunnorm_filtbetas <- getBeta(placmet_adjFunnorm_allfiltered)
placmet_adjFunnorm_filtfun_F <- placmet_adjFunnorm_filtbetas[, females$Sample_Name]
# 83 females - gives all the GSMs in either Control or PE that are female 
placmet_adjFunnorm_filtfun_M <- placmet_adjFunnorm_filtbetas[, males$Sample_Name]
# 89 males - gives all the GSMs in either Control or PE that are male 

#Calculating the average Delta Beta values between Control and PE 
Controlmetadata_all <- subset(metadata, metadata$pathology_group == "Control") 
Controlmetadata_F <- subset(females, females$pathology_group == "Control")
Controlmetadata_M <- subset(males, males$pathology_group == "Control")
# 72 control all (both males and females)
# 35 control females 
# 37 control males
PEmetadata_all <- subset(metadata, metadata$pathology_group == "PE")
PEmetadata_F <- subset(females, females$pathology_group == "PE")
PEmetadata_M <- subset(males, males$pathology_group == "PE")
# 100 PE all
# 48 PE females 
# 52 PE males

# Sex-Stratified Beta Values for PE and Control

PE_F_betas <-  as.data.frame(placmet_adjFunnorm_filtfun_F[, PEmetadata_F$Sample_Name]) #dim: 424750, 48
CONT_F_betas <- as.data.frame(placmet_adjFunnorm_filtfun_F[, Controlmetadata_F$Sample_Name]) #dim: 424750, 35

PE_M_betas <- as.data.frame(placmet_adjFunnorm_filtfun_M[, PEmetadata_M$Sample_Name])#dim: 424750, 52
CONT_M_betas <- as.data.frame(placmet_adjFunnorm_filtfun_M[, Controlmetadata_M$Sample_Name]) #dim: 424750, 37


#Cell Deconvolution - to predict placental cell composition in each sample

Cell_Decon <- function(beta){
  epidish_RPC_Group <- epidish(
    beta.m = beta[rownames(plCellCpGsThird), ],
    ref.m = plCellCpGsThird, 
    method = "RPC")

  estF_Group <- epidish_RPC_Group$estF %>% 
    as.data.frame() %>% mutate(algorithm = "RPC")
  
  return(estF_Group)}

estF_PE_F <- Cell_Decon(PE_F_betas)
estF_CONT_F <- Cell_Decon(CONT_F_betas)
estF_PE_M <- Cell_Decon(PE_M_betas)
estF_CONT_M <- Cell_Decon(CONT_M_betas)

estF_PE_F$group <- 'PE_F'
estF_CONT_F$group <- 'CONT_F'
estF_PE_M$group <- 'PE_M'
estF_CONT_M$group <- 'CONT_M'

estF <- rbind(estF_PE_F, estF_CONT_F, estF_PE_M, estF_CONT_M)




estF_diseaseandsex <- estF %>%
  mutate(pathology_group = case_when( #Creating the pathology_group variable column
    group == "PE_F" | group == "PE_M" ~ "PE",
    group == "CONT_F" | group == "CONT_M" ~ "CONT",
    .default = NA)) %>%
  mutate(sex = case_when( #Creating the sex variable column
    group == "PE_F" | group == "CONT_F" ~ "F",
    group == "PE_M" | group == "CONT_M" ~ "M",
    .default = NA))
estF_diseaseandsex$pathology_group <- as.factor(estF_diseaseandsex$pathology_group)
estF_diseaseandsex$sex <- as.factor(estF_diseaseandsex$sex)
estF_diseaseandsex$Sample_Name <- rownames(estF_diseaseandsex)
estF_dat <- merge(estF_diseaseandsex, metadata, by = "Sample_Name")
estF_dat$Fetal_Sex <- as.factor(estF_dat$Fetal_Sex)
names(estF_dat)[names(estF_dat) == "pathology_group.x"] <- "pathology_group"
cell_types <- c("Trophoblasts", "Stromal", "Hofbauer", "Endothelial", "nRBC", "Syncytiotrophoblast")

#Linear Model Comparing Cell Proportions Between Groups 
#Males: Comparing PE sample cell proportions to Control sample cell proportions
male_cell_results <- data.frame(
  Cell_Type = character(),
  Beta_PathGroup = numeric(),
  P_PathGroup = numeric(),
  stringsAsFactors = FALSE
)
male_data <- subset(estF_dat, Fetal_Sex == "M")
for(cell_type in cell_types) {
  cat("Analyzing cell types:", cell_type, "\n")
  
  #Fit linear model (comparing conception type)
  formula <- as.formula(paste(cell_type, "~ pathology_group + gestational_age + (1 | GSE_number)"))
  lmm <- lmer(formula, data = male_data)
  
  #model summary 
  model_summary <- summary(lmm)
  
  #Extract fixed effects 
  fixed_effects <-  model_summary$coefficients 
  
  #extract results for conception type 
  beta_pathology_group <- fixed_effects["pathology_groupPE", "Estimate"]
  p_pathology_group <- fixed_effects["pathology_groupPE", "Pr(>|t|)"]
  
  # Append results to the data frame
  male_cell_results <- rbind(male_cell_results, data.frame(
    Cell_Type = cell_type,
    Beta_PathGroup = beta_pathology_group,
    P_PathGroup = p_pathology_group
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
female_data <- subset(estF_dat, Fetal_Sex == "F")
for(cell_type in cell_types) {
  cat("Analyzing cell types:", cell_type, "\n")
  
  #Fit linear model (comparing conception type)
  formula <- as.formula(paste(cell_type, "~ pathology_group + gestational_age + (1 | GSE_number)"))
  lmm <- lmer(formula, data = female_data)
  
  #model summary 
  model_summary <- summary(lmm)
  
  #Extract fixed effects 
  fixed_effects <-  model_summary$coefficients 
  
  #extract results for conception type 
  beta_pathology_group <- fixed_effects["pathology_groupPE", "Estimate"]
  p_pathology_group <- fixed_effects["pathology_groupPE", "Pr(>|t|)"]
  
  # Append results to the data frame
  female_cell_results <- rbind(female_cell_results, data.frame(
    Cell_Type = cell_type,
    Beta_PathGroup = beta_pathology_group,
    P_PathGroup = p_pathology_group
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
PE_data <- subset(estF_dat, pathology_group == "PE")
for(cell_type in cell_types) {
  cat("Analyzing cell types:", cell_type, "\n")
  
  #Fit linear model (comparing conception type)
  formula <- as.formula(paste(cell_type, "~ Fetal_Sex + gestational_age + (1 | GSE_number)"))
  lmm <- lmer(formula, data = PE_data)
  
  #model summary 
  model_summary <- summary(lmm)
  
  #Extract fixed effects 
  fixed_effects <-  model_summary$coefficients 
  
  #extract results for conception type 
  beta_sex <- fixed_effects["Fetal_SexM", "Estimate"]
  p_sex <- fixed_effects["Fetal_SexM", "Pr(>|t|)"]
  
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
CONT_data <- subset(estF_dat, pathology_group == "CONT")
for(cell_type in cell_types) {
  cat("Analyzing cell types:", cell_type, "\n")
  
  #Fit linear model (comparing conception type)
  formula <- as.formula(paste(cell_type, "~ Fetal_Sex + gestational_age + (1 | GSE_number)"))
  lmm <- lmer(formula, data = CONT_data)
  
  #model summary 
  model_summary <- summary(lmm)
  
  #Extract fixed effects 
  fixed_effects <-  model_summary$coefficients 
  
  #extract results for conception type 
  beta_sex <- fixed_effects["Fetal_SexM", "Estimate"]
  p_sex <- fixed_effects["Fetal_SexM", "Pr(>|t|)"]
  
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
  formula <- as.formula(paste(cell_type, "~ Fetal_Sex + pathology_group + gestational_age + (1 | GSE_number)"))
  lmm <- lmer(formula, data = estF_dat)
  
  # Summarize the model
  model_summary <- summary(lmm)
  
  # Extract fixed effects
  fixed_effects <- model_summary$coefficients
  
  # Extract results for sex and conception type
  beta_sex <- fixed_effects["Fetal_SexM", "Estimate"]
  p_sex <- fixed_effects["Fetal_SexM", "Pr(>|t|)"]
  
  beta_pathology_group <- fixed_effects["pathology_groupPE", "Estimate"]
  p_pathology_group <- fixed_effects["pathology_groupPE", "Pr(>|t|)"]
  
  # Append results to the data frame
  cell_unstratified_results <- rbind(cell_unstratified_results, data.frame(
    Cell_Type = cell_type,
    Beta_Sex = beta_sex,
    P_Sex = p_sex,
    Beta_PathGroup = beta_pathology_group,
    P_PathGroup = p_pathology_group
  ))
}
write.csv(cell_unstratified_results, "LMM_CellType_Results.csv", row.names = FALSE)


#Plotting 
#Reshape data into long format for ggplot2
#Plot separated by study 
long_data <- estF_dat %>%
  pivot_longer(
    cols = c(Trophoblasts, Stromal, Hofbauer, Endothelial, nRBC, Syncytiotrophoblast),
    names_to = "Cell_Type",
    values_to = "Proportion"
  )
png("celltype.png", width=4000, height=1500, res=200)
ggplot(long_data, aes(x = Cell_Type, y = Proportion, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(group = group), width = 0.2, alpha = 0.1, size = 1) +
  facet_wrap(~ GSE_number, scales = "free_y") +
  labs(
    title = "Boxplot of Cell Type Proportions by Group and Sex - DNAm",
    x = "Group",
    y = "Proportion",
    fill = "Sex"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 14),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14)
  )
dev.off()
#plot with all samples in one figure
png("celltype_oneplot.png",  height = 7.5, width = 10, units = "in", res = 750)
ggplot(long_data, aes(x = Cell_Type, y = Proportion, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(group = group), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75),alpha = 0.1, size = 1) +
  labs(
    title = "Boxplot of Cell Type Proportions by Group - DNAm",
    x = "Cell Type",
    y = "Proportion",
    fill = "Group"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 14),
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    legend.position = "bottom"
  )
dev.off()
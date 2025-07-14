library(ggplot2)
library(ggsignif)
library(ggpubr)
library(tidyverse)
library(planet)
library(minfi)
library(EpiDISH)
data("plBetas")
data("plCellCpGsThird")

#Creating groups for sex stratification
metadata <- read.csv("/workspace/lab/wilsonslab/eyerk/2025_Lipid_Candidate_GSE_Info/GSE_metadata/Metadata_Sheet_lipid_preeclampsia_excluded_removed.csv")
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
placmet_adjFunnorm_allfiltered <- readRDS("/workspace/lab/wilsonslab/eyerk/2025_Thesis_Lipid_Candidate/R_entries/placmet_adjFunnorm_allfiltered.rds") #dim 329533 180
placmet_adjFunnorm_filtbetas <- getBeta(placmet_adjFunnorm_allfiltered)
placmet_adjFunnorm_filtfun_F <- placmet_adjFunnorm_filtbetas[, females$Sample_Name]
# 87 females - gives all the GSMs in either Control or PE that are female 
placmet_adjFunnorm_filtfun_M <- placmet_adjFunnorm_filtbetas[, males$Sample_Name]
# 93 males - gives all the GSMs in either Control or PE that are male 

#Calculating the average Delta Beta values between Control and PE 
Controlmetadata_all <- subset(metadata, metadata$pathology_group == "Control") 
Controlmetadata_F <- subset(females, females$pathology_group == "Control")
Controlmetadata_M <- subset(males, males$pathology_group == "Control")
# 80 control all (both males and females)
# 39 control females 
# 41 control males
PEmetadata_all <- subset(metadata, metadata$pathology_group == "PE")
PEmetadata_F <- subset(females, females$pathology_group == "PE")
PEmetadata_M <- subset(males, males$pathology_group == "PE")
# 100 PE all
# 48 PE females 
# 52 PE males

# Sex-Stratified Beta Values for PE and Control

PE_F_betas <-  as.data.frame(placmet_adjFunnorm_filtfun_F[, PEmetadata_F$Sample_Name]) #dim: 329533, 48
CONT_F_betas <- as.data.frame(placmet_adjFunnorm_filtfun_F[, Controlmetadata_F$Sample_Name]) #dim: 329533, 39

PE_M_betas <- as.data.frame(placmet_adjFunnorm_filtfun_M[, PEmetadata_M$Sample_Name])#dim: 329533, 52
CONT_M_betas <- as.data.frame(placmet_adjFunnorm_filtfun_M[, Controlmetadata_M$Sample_Name]) #dim: 329533, 41


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


# ANOVA for each cell type 
Troph.aov <- aov(Trophoblasts ~ group, data = estF)
summary(Troph.aov)
strom.aov <- aov(Stromal ~ group, data = estF)
summary(strom.aov)
hof.aov <- aov(Hofbauer ~ group, data = estF)
summary(hof.aov)
endo.aov <- aov(Endothelial ~ group, data = estF)
summary(endo.aov)
nRBC.aov <- aov(nRBC ~ group, data = estF)
summary(nRBC.aov)
syn.aov <- aov(Syncytiotrophoblast ~ group, data = estF)
summary(syn.aov)

# Bonferroni post-hoc test
troph.bonf <- pairwise.t.test(estF$Trophoblasts, estF$group, p.adjust.method = 'bonferroni')
strom.bonf <- pairwise.t.test(estF$Stromal, estF$group, p.adjust.method = 'bonferroni')
hof.bonf <- pairwise.t.test(estF$Hofbauer, estF$group, p.adjust.method = 'bonferroni')
endo.bonf <- pairwise.t.test(estF$Endothelial, estF$group, p.adjust.method = 'bonferroni')
nRBC.bonf <- pairwise.t.test(estF$nRBC, estF$group, p.adjust.method = 'bonferroni')
syn.bonf <- pairwise.t.test(estF$Syncytiotrophoblast, estF$group, p.adjust.method = 'bonferroni')

# Plotting 
cell_table <- data.frame(type = rep(c('troph', 'strom', 'hof', 'endo', 'nRBC', 'syn'), each = 4),
                         u = c(mean(estF_PE_F$Trophoblasts), mean(estF_PE_M$Trophoblasts), mean(estF_CONT_F$Trophoblasts), mean(estF_CONT_M$Trophoblasts),
                               mean(estF_PE_F$Stromal), mean(estF_PE_M$Stromal), mean(estF_CONT_F$Stromal), mean(estF_CONT_M$Stromal),
                               mean(estF_PE_F$Hofbauer), mean(estF_PE_M$Hofbauer), mean(estF_CONT_F$Hofbauer), mean(estF_CONT_M$Hofbauer),
                               mean(estF_PE_F$Endothelial), mean(estF_PE_M$Endothelial), mean(estF_CONT_F$Endothelial), mean(estF_CONT_M$Endothelial),
                               mean(estF_PE_F$nRBC), mean(estF_PE_M$nRBC), mean(estF_CONT_F$nRBC), mean(estF_CONT_M$nRBC),
                               mean(estF_PE_F$Syncytiotrophoblast), mean(estF_PE_M$Syncytiotrophoblast), mean(estF_CONT_F$Syncytiotrophoblast), mean(estF_CONT_M$Syncytiotrophoblast)),
                         s = c(sd(estF_PE_F$Trophoblasts), sd(estF_PE_M$Trophoblasts), sd(estF_CONT_F$Trophoblasts), sd(estF_CONT_M$Trophoblasts),
                               sd(estF_PE_F$Stromal), sd(estF_PE_M$Stromal), sd(estF_CONT_F$Stromal), sd(estF_CONT_M$Stromal),
                               sd(estF_PE_F$Hofbauer), sd(estF_PE_M$Hofbauer), sd(estF_CONT_F$Hofbauer), sd(estF_CONT_M$Hofbauer),
                               sd(estF_PE_F$Endothelial), sd(estF_PE_M$Endothelial), sd(estF_CONT_F$Endothelial), sd(estF_CONT_M$Endothelial),
                               sd(estF_PE_F$nRBC), sd(estF_PE_M$nRBC), sd(estF_CONT_F$nRBC), sd(estF_CONT_M$nRBC),
                               sd(estF_PE_F$Syncytiotrophoblast), sd(estF_PE_M$Syncytiotrophoblast), sd(estF_CONT_F$Syncytiotrophoblast), sd(estF_CONT_M$Syncytiotrophoblast)),
                         group = rep(c('PE_F','PE_M','CONT_F','CONT_M'), 6)
)

png(filename = "./cell_decon_anova_adjFunnorm.png", height = 7.5, width = 10, units = "in", res = 750)
ggplot(cell_table, aes(fill = group, y = u, x = type)) +
  geom_bar(position = 'dodge', stat = 'identity') +
  scale_fill_manual(values = c("#C77CFF", "#F8766D","#00BFC4", "#7CAE00")) +
  geom_errorbar(aes(ymin = u-s, ymax = u+s), width = .2, position = position_dodge(0.9)) +
  theme_classic()
dev.off()















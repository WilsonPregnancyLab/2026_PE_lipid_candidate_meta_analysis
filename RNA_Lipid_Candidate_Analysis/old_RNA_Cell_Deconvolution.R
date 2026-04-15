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
library(Biobase)
library(SingleCellExperiment)

#Bash
wget https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE182381&format=file&file=GSE182381_reference_sample.txt.gz
gunzip GSE182381_reference_sample.txt.gz

#Load RNA Read Counts
read_counts <- read.csv("/workspace/lab/wilsonslab/datalake-wilsonslab/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/read_counts/read_count_nohead.csv")
RNA_metadata <- read.csv("/workspace/lab/wilsonslab/datalake-wilsonslab/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/metadata/RNA_Lipid_Candidate_Metadata.csv")
reference_scRNAseq <- read.delim("GSE182381_reference_sample.txt")

#alter Read Counts file to be in correct format (1st column is titled 'gene names' and row values are gene name values, subsequent column names are sample names)
gtf_file <- "../genome_mapping/gencode.v48.primary_assembly.annotation.gtf"
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

# read_counts_final$RowSum <- rowSums(read_counts_final[numeric])
# dups_later <- duplicated(read_counts_final$gene_name)
# dups_earlier <- duplicated(read_counts_final$gene_name, fromLast = TRUE)
# duplicated_read_counts <- read_counts_final[dups_later | dups_earlier, ] #2111


row.names(read_counts_final_df) <- read_counts_final_df$gene_name
read_counts_final_df <- read_counts_final_df[,(2:123)]
write.csv(read_counts_final_df, "read_counts_final.csv")


# MuSiC Input Data Setup

## Reference Dataset
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

# metadata <- data.frame(labelDescription= c("sample_id", "subject_name", "cell_type"), row.names=c("sample_id", "subject_name", "cell_type"))
# sc_eset <- ExpressionSet(assayData = data.matrix(gene_exprs_mtx), phenoData =  new("AnnotatedDataFrame", data = pheno_matrix, varMetadata = metadata) )
sc_sce <- SingleCellExperiment(list(counts = gene_exprs_mtx_condensed), colData = pheno_matrix)


## Bulk Dataset
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
Est_prop_plac <-  music2_prop_t_statistics(bulk.control.mtx = bulk_control_mtx_final, bulk.case.mtx = bulk_PE_mtx_final, sc.sce = sc_sce, clusters = 'cell_id',
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

estF <- rbind(estPlac_PE_F, estPlac_CONT_F, estPlac_PE_M, estPlac_CONT_M)
est_F <- as.numeric(est_F[1:27])

# ANOVA for each cell type 
mesen.aov <- aov(Fetal.Mesenchymal.Stem.Cells ~ group, data = estF)
summary(mesen.aov)
fet.nat.kill.t.aov <- aov(Fetal.Natural.Killer.T.Cells ~ group, data = estF)
summary(fet.nat.kill.t.aov)
fet.b.aov <- aov(Fetal.B.Cells ~ group, data = estF)
summary(fet.b.aov)
fet.GZMK.nat.kill.aov <- aov(Fetal.GZMK..Natural.Killer ~ group, data = estF)
summary(fet.GZMK.nat.kill.aov)
fet.GZMB.nat.kill.aov <- aov(Fetal.GZMB..Natural.Killer ~ group, data = estF)
summary(fet.GZMB.nat.kill.aov)
fibro.aov <- aov(Fetal.Fibroblasts ~ group, data = estF)
summary(fibro.aov)
extravill.troph.aov <- aov(Fetal.Extravillous.Trophoblasts ~ group, data = estF)
summary(extravill.troph.aov)
mat.CD8.aov <- aov(Maternal.CD8..Activated.T.Cells ~ group, data = estF)
summary(mat.CD8.aov)
mat.FCGR3A.aov <- aov(Maternal.FCGR3A..Monocytes ~ group, data = estF)
summary(mat.FCGR3A.aov)
mat.CD14.troph.aov <- aov(Maternal.CD14..Monocytes ~ group, data = estF)
summary(mat.CD14.troph.aov)
mat.nat.kill.aov <- aov(Maternal.Natural.Killer.Cells ~ group, data = estF)
summary(mat.nat.kill.aov)
mat.b.aov <- aov(Maternal.B.Cells ~ group, data = estF)
summary(mat.b.aov)
mat.plas.aov <- aov(Maternal.Plasma.Cells ~ group, data = estF)
summary(mat.plas.aov)

cyto.aov <- aov(Fetal.Cytotrophoblasts ~ group, data = estF)
summary(cyto.aov)
#             Df Sum Sq Mean Sq F value   Pr(>F)
# group         3 0.1467 0.04892   12.22 5.11e-07 ***
# Residuals   118 0.4725 0.00400

p.cyto.aov <- aov(Fetal.Proliferative.Cytotrophoblasts ~ group, data = estF)
summary(p.cyto.aov)
#              Df    Sum Sq   Mean Sq F value Pr(>F)
# group         3 6.100e-06 2.033e-06   4.719 0.0038 **
# Residuals   118 5.083e-05 4.307e-07

hof.aov <- aov(Fetal.Hofbauer.Cells ~ group, data = estF)
summary(hof.aov)
#              Df Sum Sq Mean Sq F value Pr(>F)
# group         3 0.1425 0.04750   3.227 0.0251 *
# Residuals   118 1.7368 0.01472

endo.aov <- aov(Fetal.Endothelial.Cells ~ group, data = estF)
summary(endo.aov)
#              Df  Sum Sq  Mean Sq F value   Pr(>F)
# group         3 0.04244 0.014147   6.932 0.000244 ***
# Residuals   118 0.24080 0.002041

nRBC.aov <- aov(Fetal.Nucleated.Red.Blood.Cells ~ group, data = estF)
summary(nRBC.aov)
#              Df   Sum Sq   Mean Sq F value   Pr(>F)
# group         3 0.001657 0.0005523   7.242 0.000168 ***
# Residuals   118 0.008999 0.0000763

syn.aov <- aov(Fetal.Syncytiotrophoblast ~ group, data = estF)
summary(syn.aov)
#              Df Sum Sq Mean Sq F value  Pr(>F)
# group         3  1.107  0.3691   5.391 0.00164 **
# Residuals   118  8.079  0.0685

# Bonferroni post-hoc test
cyto.bonf <- pairwise.t.test(estF$Fetal.Cytotrophoblasts, estF$group, p.adjust.method = 'bonferroni')
p.cyto.bonf <- pairwise.t.test(estF$Fetal.Proliferative.Cytotrophoblasts, estF$group, p.adjust.method = 'bonferroni')
hof.bonf <- pairwise.t.test(estF$Fetal.Hofbauer.Cells, estF$group, p.adjust.method = 'bonferroni')
endo.bonf <- pairwise.t.test(estF$Fetal.Endothelial.Cells, estF$group, p.adjust.method = 'bonferroni')
nRBC.bonf <- pairwise.t.test(estF$Fetal.Nucleated.Red.Blood.Cells, estF$group, p.adjust.method = 'bonferroni')
syn.bonf <- pairwise.t.test(estF$Fetal.Syncytiotrophoblast, estF$group, p.adjust.method = 'bonferroni')

# Plotting 
cell_table <- data.frame(type = rep(c('cyto', 'p.cyto', 'hof', 'endo', 'nRBC', 'syn'), each = 4),
                         u = c(mean(estF_PE_F$Fetal.Cytotrophoblasts), mean(estF_PE_M$Fetal.Cytotrophoblasts), mean(estF_CONT_F$Fetal.Cytotrophoblasts), mean(estF_CONT_M$Fetal.Cytotrophoblasts),
                               mean(estF_PE_F$Fetal.Proliferative.Cytotrophoblasts), mean(estF_PE_M$Fetal.Proliferative.Cytotrophoblasts), mean(estF_CONT_F$Fetal.Proliferative.Cytotrophoblasts), mean(estF_CONT_M$Fetal.Proliferative.Cytotrophoblasts),
                               mean(estF_PE_F$Fetal.Hofbauer.Cells), mean(estF_PE_M$Fetal.Hofbauer.Cells), mean(estF_CONT_F$Fetal.Hofbauer.Cells), mean(estF_CONT_M$Fetal.Hofbauer.Cells),
                               mean(estF_PE_F$Fetal.Endothelial.Cells), mean(estF_PE_M$Fetal.Endothelial.Cells), mean(estF_CONT_F$Fetal.Endothelial.Cells), mean(estF_CONT_M$Fetal.Endothelial.Cells),
                               mean(estF_PE_F$Fetal.Nucleated.Red.Blood.Cells), mean(estF_PE_M$Fetal.Nucleated.Red.Blood.Cells), mean(estF_CONT_F$Fetal.Nucleated.Red.Blood.Cells), mean(estF_CONT_M$Fetal.Nucleated.Red.Blood.Cells),
                               mean(estF_PE_F$Fetal.Syncytiotrophoblast), mean(estF_PE_M$Fetal.Syncytiotrophoblast), mean(estF_CONT_F$Fetal.Syncytiotrophoblast), mean(estF_CONT_M$Fetal.Syncytiotrophoblast)),
                         s = c(sd(estF_PE_F$Fetal.Cytotrophoblasts), sd(estF_PE_M$Fetal.Cytotrophoblasts), sd(estF_CONT_F$Fetal.Cytotrophoblasts), sd(estF_CONT_M$Fetal.Cytotrophoblasts),
                               sd(estF_PE_F$Fetal.Proliferative.Cytotrophoblasts), sd(estF_PE_M$Fetal.Proliferative.Cytotrophoblasts), sd(estF_CONT_F$Fetal.Proliferative.Cytotrophoblasts), sd(estF_CONT_M$Fetal.Proliferative.Cytotrophoblasts),
                               sd(estF_PE_F$Fetal.Hofbauer.Cells), sd(estF_PE_M$Fetal.Hofbauer.Cells), sd(estF_CONT_F$Fetal.Hofbauer.Cells), sd(estF_CONT_M$Fetal.Hofbauer.Cells),
                               sd(estF_PE_F$Fetal.Endothelial.Cells), sd(estF_PE_M$Fetal.Endothelial.Cells), sd(estF_CONT_F$Fetal.Endothelial.Cells), sd(estF_CONT_M$Fetal.Endothelial.Cells),
                               sd(estF_PE_F$Fetal.Nucleated.Red.Blood.Cells), sd(estF_PE_M$Fetal.Nucleated.Red.Blood.Cells), sd(estF_CONT_F$Fetal.Nucleated.Red.Blood.Cells), sd(estF_CONT_M$Fetal.Nucleated.Red.Blood.Cells),
                               sd(estF_PE_F$Fetal.Syncytiotrophoblast), sd(estF_PE_M$Fetal.Syncytiotrophoblast), sd(estF_CONT_F$Fetal.Syncytiotrophoblast), sd(estF_CONT_M$Fetal.Syncytiotrophoblast)),
                         group = rep(c('PE_F','PE_M','CONT_F','CONT_M'), 6)
)

png(filename = "./RNA_cell_decon_anova.png", height = 7.5, width = 10, units = "in", res = 750)
ggplot(cell_table, aes(fill = group, y = u, x = type)) +
  geom_bar(position = 'dodge', stat = 'identity') +
  scale_fill_manual(values = c("#C77CFF", "#F8766D","#00BFC4", "#7CAE00")) +
  geom_errorbar(aes(ymin = u-s, ymax = u+s), width = .2, position = position_dodge(0.9)) +
  theme_classic()
dev.off()


cols <-c("cyto" = "cadetblue2", "p.cyto" = "lightsalmon1", "hof" = "palegreen2", "endo" = "goldenrod1",
          "nRBC"="steelblue3", "syn" = "plum2")
          
png(filename = "./RNA_cell_decon_anova_dotty.png", height = 7.5, width = 10, units = "in", res = 750)

ggplot(cell_table, aes(x=type, y=u, fill=type)) + xlab('')+
  geom_jitter(width=0.25,alpha=0.8)+ylab('Cell Type Proportions')+theme_bw()+
  stat_summary(fun = median,
               geom = "crossbar", width = 0.5,linewidth=0.5,color='gray36')+
  facet_grid(.~group)+
  theme(plot.title = element_text(hjust = 0.5, size=12),
        axis.text.x = element_text(size=12,angle = 45,hjust=1),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        axis.line = element_line(colour = "black"),
        strip.text.x = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')+
  scale_fill_manual(values=cols)
dev.off()


# plot estimated cell type proportions
prop_all <- Final_est_plac_prop %>%
            tibble::rownames_to_column(var = "sampleID") %>%
            pivot_longer(
    cols = -sampleID, 
    names_to = "celltype",          
    values_to = "proportion"       
  ) %>%
  mutate(celltype = factor(celltype),
    proportion = as.numeric(proportion)
  )

prop_all <- as.data.frame(prop_all)
prop_all$group <- ifelse(prop_all$sampleID %in% rownames(metadata_keep_eset_control), yes = "Control", no = "PE")
cells_to_keep <- c("Fetal.Cytotrophoblasts", "Fetal.Proliferative.Cytotrophoblasts", "Fetal.Hofbauer.Cells","Fetal.Endothelial.Cells","Fetal.Nucleated.Red.Blood.Cells","Fetal.Syncytiotrophoblast")
prop_select <- prop_all[prop_all$celltype %in% cells_to_keep,]

cols <-c("Fetal.Cytotrophoblasts" = "cadetblue2", "Fetal.Proliferative.Cytotrophoblasts" = "lightsalmon1", "Fetal.Hofbauer.Cells" = "palegreen2", "Fetal.Endothelial.Cells" = "goldenrod1",
          "Fetal.Nucleated.Red.Blood.Cells"="steelblue3", "Fetal.Syncytiotrophoblast" = "plum2")

png(filename = "./RNA_cell_decon_anova_dotty.png", height = 7.5, width = 10, units = "in", res = 750)

ggplot(prop_select, aes(x=celltype, y=proportion, color=celltype)) + xlab('')+
  geom_jitter(width=0.25,alpha=0.8)+ylab('Cell Type Proportions')+theme_bw()+
  stat_summary(fun = median,
               geom = "crossbar", width = 0.5,size=0.5,color='gray36')+
  facet_grid(.~group)+
  theme(plot.title = element_text(hjust = 0.5, size=12),
        axis.text.x = element_text(size=12,angle = 45,hjust=1),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        axis.line = element_line(colour = "black"),
        strip.text.x = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')+
  scale_color_manual(values=cols)

dev.off()


cyto.aov <- pairwise.t.test(estF$Fetal.Cytotrophoblasts, estF$group, p.adjust.method = 'bonferroni')
p.cyto.aov <- pairwise.t.test(estF$Fetal.Proliferative.Cytotrophoblasts, estF$group, p.adjust.method = 'bonferroni')
hof.bonf <- pairwise.t.test(estF$Fetal.Hofbauer.Cells, estF$group, p.adjust.method = 'bonferroni')
endo.bonf <- pairwise.t.test(estF$Fetal.Endothelial.Cells, estF$group, p.adjust.method = 'bonferroni')
nRBC.bonf <- pairwise.t.test(estF$Fetal.Nucleated.Red.Blood.Cells, estF$group, p.adjust.method = 'bonferroni')
syn.bonf <- pairwise.t.test(estF$Fetal.Syncytiotrophoblast, estF$group, p.adjust.method = 'bonferroni')











"Fetal.Mesenchymal.Stem.Cells", 
#"Fetal.CD14..Monocytes", 
#"Fetal.CD8..Activated.T.Cells", 
#"Fetal.Naive.CD4..T.Cells", 
#"Fetal.Naive.CD8..T.Cells", 
"Fetal.Natural.Killer.T.Cells", 
"Fetal.B.Cells", 
"Fetal.GZMK..Natural.Killer",
#"Fetal.Memory.CD4..T.Cells",  
"Fetal.Hofbauer.Cells Fetal",
#"Fetal.Plasmacytoid.Dendritic.Cells", 
"Fetal.GZMB..Natural.Killer", 
"Fetal.Endothelial.Cells",
"Fetal.Syncytiotrophoblast",
"Fetal.Fibroblasts", 
"Fetal.Cytotrophoblasts", 
"Fetal.Proliferative.Cytotrophoblasts",  
"Fetal.Nucleated.Red.Blood.Cells",
"Maternal.CD8..Activated.T.Cells", 
#"Maternal.Naive.CD4..T.Cells", 
"Maternal.FCGR3A..Monocytes", 
"Maternal.CD14..Monocytes", 
"Maternal.Natural.Killer.Cells", 
"Maternal.B.Cells", 
"Maternal.Plasma.Cells", 
#"Maternal.Naive.CD8..T.Cells", 
"Fetal.Extravillous.Trophoblasts"












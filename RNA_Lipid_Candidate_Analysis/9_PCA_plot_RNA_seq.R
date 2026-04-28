library(DESeq2)
library(limma)
library(ggplot2) #version 3.5.1
library(gridExtra) #version 2.3
library(ggrepel) 

#Before Controlling for Batch Effects 
#Load DESeqDataSet object from RNA DESeq Analysis (rerun first part of "8a_DESeq2_Lip_Candidate_differential_exp_analysis.R")
#Colour each PCA by group, fetal-sex, and disease group

filt_read_counts_autosomes <- filt_read_counts[rownames(filt_read_counts) %in% autosomal_lipids$gene_id,]
filt_read_counts_chrX <- filt_read_counts [rownames(filt_read_counts) %in% chrX_lipids$gene_id,]
filt_read_counts_chrY <- filt_read_counts [rownames(filt_read_counts) %in% chrY_lipids$gene_id,]

#Rownames in metadata match with colnames of filt_read_counts_ComBat, and in same order
rownames(RNA_Metadata_ex_rem) <- RNA_Metadata_ex_rem$Run
all(colnames(filt_read_counts_ComBat) %in% rownames(RNA_Metadata_ex_rem)) #TRUE
all(colnames(filt_read_counts_ComBat) == rownames(RNA_Metadata_ex_rem)) #TRUE

#Step 3: Construct DESeqDataSet object, removed "instrument" from design because confound perfectly with GSE_number
RNA_WG_DESeq_auto_ucor <- DESeqDataSetFromMatrix(countData = filt_read_counts_autosomes, colData = RNA_Metadata_ex_rem, design = ~ disease_group + GSE_number + predicted_fetal_sex)

#Plot PCA
vsd_Lip_RNA_batch_uncorrected <- vst(RNA_WG_DESeq_auto_ucor, blind = FALSE)

PCA_batch_uncorrected_study <- plotPCA(vsd_Lip_RNA_batch_uncorrected, intgroup = "GSE_number")
PCA_batch_uncorrected_study_plot <- PCA_batch_uncorrected_study + 
    labs(title = "PCA of Batch Uncorrected RNA-Seq Data - Study") +
  theme_bw() +
  theme(plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 14)) +  
    scale_y_continuous(breaks = seq(-25, 25, by = 10), limits = c(-25, 25)) +
    scale_x_continuous(breaks = seq(-30, 30, by = 10), limits = c(-30, 30)) +
    stat_ellipse(geom = "polygon",
               aes(fill = group), 
               alpha = 0.25)
png("./PCA_batch_uncorrected_study.png", height = 10, width = 18, units = "in", res = 300)
print(PCA_batch_uncorrected_study_plot, nrow = 1)
dev.off()

PCA_batch_uncorrected_fetal_sex <- plotPCA(vsd_Lip_RNA_batch_uncorrected, intgroup = "predicted_fetal_sex")
PCA_batch_uncorrected_fetal_sex_plot <- PCA_batch_uncorrected_fetal_sex + 
    labs(title = "PCA of Batch Uncorrected RNA-Seq Data - Fetal Sex") +
  theme_bw() +
  theme(plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 14)) +
    scale_y_continuous(breaks = seq(-20, 20, by = 10), limits = c(-20, 20)) +
    scale_x_continuous(breaks = seq(-50, 50, by = 10), limits = c(-50, 50)) +
    stat_ellipse(geom = "polygon",
               aes(fill = group), 
               alpha = 0.25)
png("./PCA_batch_uncorrected_fetalsex.png", height = 10, width = 18, units = "in", res = 300)
grid.arrange(PCA_batch_uncorrected_fetal_sex_plot, nrow = 1)
dev.off()

PCA_batch_uncorrected_disease_group <- plotPCA(vsd_Lip_RNA_batch_uncorrected, intgroup = "disease_group")
PCA_batch_uncorrected_disease_group_plot <- PCA_batch_uncorrected_disease_group + 
    labs(title = "PCA of Batch Uncorrected RNA-Seq Data - Disease Group") +
  theme_bw() +
  theme(plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 14)) +
    scale_y_continuous(breaks = seq(-20, 20, by = 10), limits = c(-20, 20)) +
    scale_x_continuous(breaks = seq(-50, 50, by = 10), limits = c(-50, 50)) +
    stat_ellipse(geom = "polygon",
               aes(fill = group), 
               alpha = 0.25)

png("./PCA_batch_uncorrected_disease_group.png", height = 10, width = 18, units = "in", res = 300)
grid.arrange(PCA_batch_uncorrected_disease_group_plot, nrow = 1)
dev.off()


#After Controlling for Batch Effects 
vsd_Lip_RNA_batch_corrected <- vst(RNA_WG_DESeq_auto, blind = FALSE)
SV_numbers <- c("SV1", "SV2", "SV3", "SV4", "SV5", "SV6")
sv_matrix <- as.matrix(metadata_combined [, colnames(metadata_combined) %in% SV_numbers])

batch_corrected_matrix <- limma::removeBatchEffect(assay(vsd_Lip_RNA_batch_corrected), covariates=sv_matrix, design = model.matrix(~disease_group + predicted_fetal_sex, data = metadata_combined))
assay(vsd_Lip_RNA_batch_corrected) <- batch_corrected_matrix

PCA_batch_corrected_study <- plotPCA(vsd_Lip_RNA_batch_corrected, intgroup="GSE_number", returnData = TRUE)
PCA_batch_corrected_study_plot <- PCA_batch_corrected_study + 
    labs(title = "PCA of Batch Corrected RNA-Seq Data - Study") +
  theme_bw() +
  theme(plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 14)) +
        scale_y_continuous(breaks = seq(-25, 25, by = 10), limits = c(-25, 25)) +
    scale_x_continuous(breaks = seq(-30, 30, by = 10), limits = c(-30, 30)) +
    stat_ellipse(geom = "polygon",
               aes(fill = group), 
               alpha = 0.25)

png("./PCA_batch_corrected_study.png", height = 10, width = 18, units = "in", res = 300)
grid.arrange(PCA_batch_corrected_study_plot, nrow = 1)
dev.off()

PCA_batch_corrected_fetal_sex <- plotPCA(vsd_Lip_RNA_batch_corrected, intgroup = "predicted_fetal_sex")
PCA_batch_corrected_fetal_sex_plot <- PCA_batch_corrected_fetal_sex + 
    labs(title = "PCA of Batch Corrected RNA-Seq Data - Fetal Sex") +
  theme_bw() +
  theme(plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 14)) +
    scale_y_continuous(breaks = seq(-20, 20, by = 10), limits = c(-20, 20)) +
    scale_x_continuous(breaks = seq(-50, 50, by = 10), limits = c(-50, 50)) +
    stat_ellipse(geom = "polygon",
               aes(fill = group), 
               alpha = 0.25)

png("./PCA_batch_corrected_fetalsex.png", height = 10, width = 18, units = "in", res = 300)
grid.arrange(PCA_batch_corrected_fetal_sex_plot, nrow = 1)
dev.off()

PCA_batch_corrected_disease_group <- plotPCA(vsd_Lip_RNA_batch_corrected, intgroup = "disease_group")
PCA_batch_corrected_disease_group_plot <- PCA_batch_corrected_disease_group + 
    labs(title = "PCA of Batch Corrected RNA-Seq Data - Disease Group") +
  theme_bw() +
  theme(plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 14)) +
    scale_y_continuous(breaks = seq(-20, 20, by = 10), limits = c(-20, 20)) +
    scale_x_continuous(breaks = seq(-50, 50, by = 10), limits = c(-50, 50)) +
    stat_ellipse(geom = "polygon",
               aes(fill = group), 
               alpha = 0.25)
png("./PCA_batch_corrected_disease_group.png", height = 10, width = 18, units = "in", res = 300)
grid.arrange(PCA_batch_corrected_disease_group_plot, nrow = 1)
dev.off()




png("./PCA_batch_study.png", height = 8, width = 20, units = "in", res = 300)
grid.arrange(PCA_batch_uncorrected_study_plot, PCA_batch_corrected_study_plot, nrow = 1)
dev.off()

png("./PCA_batch_fetal_sex.png", height = 8, width = 20, units = "in", res = 300)
grid.arrange(PCA_batch_uncorrected_fetal_sex_plot, PCA_batch_corrected_fetal_sex_plot, nrow = 1)
dev.off()

png("./PCA_batch_disease_group.png", height = 8, width = 20, units = "in", res = 300)
grid.arrange(PCA_batch_uncorrected_disease_group_plot, PCA_batch_corrected_disease_group_plot, nrow = 1)
dev.off()

#Removed GSE234729 and redid differential expression analysis because samples clustered into 2 separate clusters and didn't recluster after normalization

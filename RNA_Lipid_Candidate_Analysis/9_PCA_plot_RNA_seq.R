library(DESeq2)
library(limma)
library(ggplot2) #version 3.5.1
library(gridExtra) #version 2.3
library(ggrepel) 

#Before Controlling for Batch Effects 
#Load DESeqDataSet object from RNA DESeq Analysis (rerun first part of "8a_DESeq2_Lip_Candidate_differential_exp_analysis.R")
#Colour each PCA by group, fetal-sex, and disease group
vsd_Lip_RNA_batch_uncorrected <- vst(RNA_Lip_DESeq_auto, blind = FALSE)

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
grid.arrange(PCA_batch_uncorrected_study_plot, nrow = 1)
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

vst_Lip_RNA_batch_uncorrected_matrix <- assay(vsd_Lip_RNA_batch_uncorrected)

vst_Lip_RNA_batch_corrected <- limma::removeBatchEffect(vst_Lip_RNA_batch_uncorrected_matrix, batch=vsd_Lip_RNA_batch_uncorrected$GSE_number, batch2=vsd_Lip_RNA_batch_uncorrected$predicted_fetal_sex, design = model.matrix(~ vsd_Lip_RNA_batch_uncorrected$disease_group))

vsd_Lip_RNA_batch_corrected <- vsd_Lip_RNA_batch_uncorrected

assay(vsd_Lip_RNA_batch_corrected) <- vst_Lip_RNA_batch_corrected

PCA_batch_corrected_study <- plotPCA(vsd_Lip_RNA_batch_corrected, intgroup="GSE_number")
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

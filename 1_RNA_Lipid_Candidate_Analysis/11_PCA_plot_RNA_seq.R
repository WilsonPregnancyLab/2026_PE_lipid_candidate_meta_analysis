library(DESeq2)
library(limma)
library(ggplot2) #version 3.5.1
library(gridExtra) #version 2.3
library(ggrepel) 

#Before Controlling for Batch Effects (n = 205; GSE143953, GSE148241, GSE186257, GSE234729)
#Load DESeqDataSet object from RNA DESeq Analysis (rerun first part of "8a_DESeq2_Lip_Candidate_differential_exp_analysis.R")
#Colour each PCA by group, fetal-sex, and disease group


RNA_metadata <- read.csv("../metadata/RNA_Lipid_Candidate_Metadata.csv")

metadata_4studies <- RNA_metadata[RNA_metadata$GSE_number %in% c("GSE143953","GSE148241","GSE186257","GSE234729"),]
write.csv(metadata_4studies, "metadata_4studies.csv")

#unexclude only GSE234729
metadata_4studies <- read.csv("metadata_4studies.csv")

#Condense read_counts to only exclude "exclude" samples from metadata sheet, and GSE234729, and separate by chromosome
keep_GSE <- metadata_4studies[!metadata_4studies$exclude == "exclude", ]
filt_read_counts_4studies <- read_counts[, colnames(read_counts) %in% keep_GSE$Run] #78894 202
keep_GSE <- keep_GSE[keep_GSE$Run %in% colnames(read_counts),] #202, 25

filt_read_counts_4studies_autosomes <- filt_read_counts_4studies[rownames(filt_read_counts_4studies) %in% autosomal_lipids$gene_id,]
filt_read_counts_4studies_chrX <- filt_read_counts_4studies [rownames(filt_read_counts_4studies) %in% chrX_lipids$gene_id,]
filt_read_counts_4studies_chrY <- filt_read_counts_4studies [rownames(filt_read_counts_4studies) %in% chrY_lipids$gene_id,]

#Rownames in metadata match with colnames of filt_read_counts_4studies, and in same order
rownames(keep_GSE) <- keep_GSE$Run
all(colnames(filt_read_counts_4studies) %in% rownames(keep_GSE)) #TRUE
all(colnames(filt_read_counts_4studies) == rownames(keep_GSE)) #TRUE

#Step 3: Construct DESeqDataSet object, removed "instrument" from design because confound perfectly with GSE_number
RNA_WG_DESeq_auto_ucor_4studies <- DESeqDataSetFromMatrix(countData = filt_read_counts_4studies_autosomes, colData = keep_GSE, design = ~ disease_group + GSE_number + predicted_fetal_sex)

#Plot PCA
vsd_Lip_RNA_batch_uncorrected_4studies <- vst(RNA_WG_DESeq_auto_ucor_4studies, blind = FALSE)

PCA_batch_uncorrected_study_4studies <- plotPCA(vsd_Lip_RNA_batch_uncorrected_4studies, intgroup = "GSE_number")
PCA_batch_uncorrected_study_plot_4studies <- PCA_batch_uncorrected_study_4studies + 
    labs(title = "PCA of Batch Uncorrected RNA-Seq Data - Study") +
  theme_classic() +
  theme(plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 14),  
        panel.grid = element_blank()) + 
    scale_y_continuous(breaks = seq(-25, 25, by = 10), limits = c(-25, 25)) +
    scale_x_continuous(breaks = seq(-30, 30, by = 10), limits = c(-30, 30)) +
    geom_vline(xintercept = c(-0.05,0.05), col = "black", linetype = "dashed", linewidth = 0.75) +
    geom_hline(yintercept = c(-log10(0.05)), col = "black", linetype = "dashed", linewidth = 0.75) +
    #stat_ellipse(geom = "polygon",
     #          aes(fill = group), 
      #         alpha = 0.25)
png("./PCA_batch_uncorrected_study_4studies.png", height = 10, width = 18, units = "in", res = 300)
print(PCA_batch_uncorrected_study_plot_4studies, nrow = 1)
dev.off()

#batch corrected
vsd_Lip_RNA_batch_corrected_4studies <- vst(RNA_WG_DESeq_auto_ucor_4studies, blind = FALSE)
# SV_numbers <- c("SV1", "SV2", "SV3", "SV4", "SV5", "SV6")
# sv_matrix <- as.matrix(RNA_Metadata_ex_rem [, colnames(RNA_Metadata_ex_rem) %in% SV_numbers])

batch_corrected_matrix_4studies <- limma::removeBatchEffect(assay(vsd_Lip_RNA_batch_corrected_4studies),  batch = keep_GSE$GSE_number, design = model.matrix(~disease_group + predicted_fetal_sex, data = keep_GSE))
assay(vsd_Lip_RNA_batch_corrected_4studies) <- batch_corrected_matrix_4studies

PCA_batch_corrected_study_4studies <- plotPCA(vsd_Lip_RNA_batch_corrected_4studies, intgroup="GSE_number")
PCA_batch_corrected_study_plot_4studies <- PCA_batch_corrected_study_4studies + 
    labs(title = "PCA of Batch Corrected RNA-Seq Data - Study") +
  theme_classic() +
  theme(plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 14),
         panel.grid = element_blank()) +
        scale_y_continuous(breaks = seq(-25, 25, by = 10), limits = c(-25, 25)) +
    scale_x_continuous(breaks = seq(-30, 30, by = 10), limits = c(-30, 30)) +
        geom_vline(xintercept = c(-0.05,0.05), col = "black", linetype = "dashed", linewidth = 0.75) +
    geom_hline(yintercept = c(-log10(0.05)), col = "black", linetype = "dashed", linewidth = 0.75) +
    # stat_ellipse(geom = "polygon",
    #            aes(fill = group), 
    #            alpha = 0.25)

png("./PCA_batch_corrected_study_4studies.png", height = 10, width = 18, units = "in", res = 300)
grid.arrange(PCA_batch_corrected_study_plot_4studies, nrow = 1)
dev.off()

png("./PCA_batch_study_4studies.png", height = 8, width = 20, units = "in", res = 300)
grid.arrange(PCA_batch_uncorrected_study_plot_4studies, PCA_batch_corrected_study_plot_4studies, nrow = 1)
dev.off()


#Before Controlling for Batch Effects (n = 92; GSE143953, GSE148241, GSE186257)
#Load DESeqDataSet object from RNA DESeq Analysis (rerun first part of "8a_DESeq2_Lip_Candidate_differential_exp_analysis.R")
#Colour each PCA by group, fetal-sex, and disease group

filt_read_counts_autosomes <- filt_read_counts[rownames(filt_read_counts) %in% autosomal_lipids$gene_id,]
filt_read_counts_chrX <- filt_read_counts [rownames(filt_read_counts) %in% chrX_lipids$gene_id,]
filt_read_counts_chrY <- filt_read_counts [rownames(filt_read_counts) %in% chrY_lipids$gene_id,]

#Rownames in metadata match with colnames of filt_read_counts_ComBat, and in same order
rownames(RNA_Metadata_ex_rem) <- RNA_Metadata_ex_rem$Run
all(colnames(filt_read_counts) %in% rownames(RNA_Metadata_ex_rem)) #TRUE
all(colnames(filt_read_counts) == rownames(RNA_Metadata_ex_rem)) #TRUE

#Step 3: Construct DESeqDataSet object, removed "instrument" from design because confound perfectly with GSE_number
RNA_WG_DESeq_auto_ucor <- DESeqDataSetFromMatrix(countData = filt_read_counts_autosomes, colData = RNA_Metadata_ex_rem, design = ~ disease_group + GSE_number + predicted_fetal_sex)

#Plot PCA
vsd_Lip_RNA_batch_uncorrected <- vst(RNA_WG_DESeq_auto_ucor, blind = FALSE)

PCA_batch_uncorrected_study <- plotPCA(vsd_Lip_RNA_batch_uncorrected, intgroup = "GSE_number")
PCA_batch_uncorrected_study_plot <- PCA_batch_uncorrected_study + 
    labs(title = "PCA of Batch Uncorrected RNA-Seq Data - Study") +
  theme_classic() +
  theme(plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 14),  
        panel.grid = element_blank()) + 
    scale_y_continuous(breaks = seq(-25, 25, by = 10), limits = c(-25, 25)) +
    scale_x_continuous(breaks = seq(-30, 30, by = 10), limits = c(-30, 30)) +
    geom_vline(xintercept = c(-0.05,0.05), col = "black", linetype = "dashed", linewidth = 0.75) +
    geom_hline(yintercept = c(-log10(0.05)), col = "black", linetype = "dashed", linewidth = 0.75) +
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
        legend.text = element_text(size = 14),
         panel.grid = element_blank()) +
    scale_y_continuous(breaks = seq(-20, 20, by = 10), limits = c(-20, 20)) +
    scale_x_continuous(breaks = seq(-50, 50, by = 10), limits = c(-50, 50)) +
    geom_vline(xintercept = c(-0.05,0.05), col = "black", linetype = "dashed", linewidth = 0.75) +
    geom_hline(yintercept = c(-log10(0.05)), col = "black", linetype = "dashed", linewidth = 0.75) +
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
        legend.text = element_text(size = 14),
         panel.grid = element_blank()) +
    scale_y_continuous(breaks = seq(-20, 20, by = 10), limits = c(-20, 20)) +
    scale_x_continuous(breaks = seq(-50, 50, by = 10), limits = c(-50, 50)) +
    geom_vline(xintercept = c(-0.05,0.05), col = "black", linetype = "dashed", linewidth = 0.75) +
    geom_hline(yintercept = c(-log10(0.05)), col = "black", linetype = "dashed", linewidth = 0.75) +
    stat_ellipse(geom = "polygon",
               aes(fill = group), 
               alpha = 0.25)

png("./PCA_batch_uncorrected_disease_group.png", height = 10, width = 18, units = "in", res = 300)
grid.arrange(PCA_batch_uncorrected_disease_group_plot, nrow = 1)
dev.off()


#After Controlling for Batch Effects 
vsd_Lip_RNA_batch_corrected <- vst(RNA_WG_DESeq_auto, blind = FALSE)
# SV_numbers <- c("SV1", "SV2", "SV3", "SV4", "SV5", "SV6")
# sv_matrix <- as.matrix(RNA_Metadata_ex_rem [, colnames(RNA_Metadata_ex_rem) %in% SV_numbers])

batch_corrected_matrix <- limma::removeBatchEffect(assay(vsd_Lip_RNA_batch_corrected),  batch = RNA_Metadata_ex_rem$GSE_number, design = model.matrix(~disease_group + predicted_fetal_sex, data = RNA_Metadata_ex_rem))
assay(vsd_Lip_RNA_batch_corrected) <- batch_corrected_matrix

PCA_batch_corrected_study <- plotPCA(vsd_Lip_RNA_batch_corrected, intgroup="GSE_number")
PCA_batch_corrected_study_plot <- PCA_batch_corrected_study + 
    labs(title = "PCA of Batch Corrected RNA-Seq Data - Study") +
  theme_classic() +
  theme(plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 14),
         panel.grid = element_blank()) +
        scale_y_continuous(breaks = seq(-25, 25, by = 10), limits = c(-25, 25)) +
    scale_x_continuous(breaks = seq(-30, 30, by = 10), limits = c(-30, 30)) +
        geom_vline(xintercept = c(-0.05,0.05), col = "black", linetype = "dashed", linewidth = 0.75) +
    geom_hline(yintercept = c(-log10(0.05)), col = "black", linetype = "dashed", linewidth = 0.75) +
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
        legend.text = element_text(size = 14),
         panel.grid = element_blank()) +
    scale_y_continuous(breaks = seq(-20, 20, by = 10), limits = c(-20, 20)) +
    scale_x_continuous(breaks = seq(-50, 50, by = 10), limits = c(-50, 50)) +
    geom_vline(xintercept = c(-0.05,0.05), col = "black", linetype = "dashed", linewidth = 0.75) +
    geom_hline(yintercept = c(-log10(0.05)), col = "black", linetype = "dashed", linewidth = 0.75) +
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
        legend.text = element_text(size = 14),
         panel.grid = element_blank()) +
    scale_y_continuous(breaks = seq(-20, 20, by = 10), limits = c(-20, 20)) +
    scale_x_continuous(breaks = seq(-50, 50, by = 10), limits = c(-50, 50)) +
    geom_vline(xintercept = c(-0.05,0.05), col = "black", linetype = "dashed", linewidth = 0.75) +
    geom_hline(yintercept = c(-log10(0.05)), col = "black", linetype = "dashed", linewidth = 0.75) +
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

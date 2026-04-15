
# ##Saw some weirdness with GSE234729, going to do PCA just of normalized GSE234729 and the rest of the 4 GSEs (GSE143953, GSE148241, GSE186257, GSE218039)
# GSE234729_Samples <- RNA_Metadata_ex_rem[RNA_Metadata_ex_rem$GSE_number == "GSE234729",]

# filt_read_counts_GSE234729 <- filt_read_counts_autosomes[, (colnames(filt_read_counts_autosomes) %in% rownames(GSE234729_Samples))]

# all(colnames(filt_read_counts_GSE234729) %in% rownames(GSE234729_Samples)) #TRUE
# all(colnames(filt_read_counts_GSE234729) == rownames(GSE234729_Samples)) #TRUE

# RNA_WG_DESeq_GSE234729 <- DESeqDataSetFromMatrix(countData = filt_read_counts_GSE234729, colData = GSE234729_Samples, design = ~ disease_group + predicted_fetal_sex)

# pre_filter <- rowSums(counts(RNA_WG_DESeq_GSE234729)) >=10
# RNA_WG_DESeq_GSE234729 <- RNA_WG_DESeq_GSE234729[pre_filter,]
# RNA_Lip_DESeq_GSE234729 <- RNA_WG_DESeq_GSE234729[rownames(RNA_WG_DESeq_GSE234729) %in% lipid_ensembl_list$gene_id]

# RNA_Lip_DESeq_GSE234729$disease_group <- relevel(RNA_Lip_DESeq_GSE234729$disease_group, ref = "Control")
# RNA_Lip_DESeq_GSE234729 <- DESeq(RNA_Lip_DESeq_GSE234729)

# vsd_Lip_RNA_batch_uncorrected_GSE234729 <- vst(RNA_Lip_DESeq_GSE234729, blind = FALSE)


# PCA_batch_uncorrected_GSE234729_disease_group <- plotPCA(vsd_Lip_RNA_batch_uncorrected_GSE234729, intgroup="disease_group")

# png("./PCA_batch_uncorrected_GSE234729_disease_group.png", height = 9, width = 15, units = "in", res = 300)
# grid.arrange(PCA_batch_uncorrected_GSE234729_disease_group, nrow = 1)
# dev.off()

# PCA_batch_uncorrected_GSE234729_ethnicity <- plotPCA(vsd_Lip_RNA_batch_uncorrected_GSE234729, intgroup="maternal_ethnicity")

# png("./PCA_batch_uncorrected_GSE234729_ethnicity.png", height = 9, width = 15, units = "in", res = 300)
# grid.arrange(PCA_batch_uncorrected_GSE234729_ethnicity, nrow = 1)
# dev.off()


# #GSE234729 batch corrected
# vst_Lip_RNA_batch_uncorrected_matrix_GSE234729 <- assay(vsd_Lip_RNA_batch_uncorrected_GSE234729)

# vst_Lip_RNA_batch_corrected_GSE234729 <- limma::removeBatchEffect(vst_Lip_RNA_batch_uncorrected_matrix_GSE234729, batch=vsd_Lip_RNA_batch_uncorrected_GSE234729$predicted_fetal_sex, design = model.matrix(~ vsd_Lip_RNA_batch_uncorrected_GSE234729$disease_group))

# vsd_Lip_RNA_batch_corrected_GSE234729 <- vsd_Lip_RNA_batch_uncorrected_GSE234729

# assay(vsd_Lip_RNA_batch_corrected_GSE234729) <- vst_Lip_RNA_batch_corrected_GSE234729

# PCA_batch_corrected_GSE234729_disease_group <- plotPCA(vsd_Lip_RNA_batch_corrected_GSE234729, intgroup="disease_group")

# png("./PCA_batch_corrected_GSE234729_disease_group.png", height = 9, width = 15, units = "in", res = 300)
# grid.arrange(PCA_batch_corrected_GSE234729_disease_group, nrow = 1)
# dev.off()


# #Rest of 4 GSEs
# Remaining_GSE <- c("GSE143953", "GSE148241", "GSE186257", "GSE218039")
# Non_GSE234729_Samples <- RNA_Metadata_ex_rem[RNA_Metadata_ex_rem$GSE_number %in% Remaining_GSE,]

# filt_read_counts_4 <- filt_read_counts_autosomes[, (colnames(filt_read_counts_autosomes) %in% rownames(Non_GSE234729_Samples))]

# all(colnames(filt_read_counts_4) %in% rownames(Non_GSE234729_Samples)) #TRUE
# all(colnames(filt_read_counts_4) == rownames(Non_GSE234729_Samples)) #TRUE

# RNA_WG_DESeq_4 <- DESeqDataSetFromMatrix(countData = filt_read_counts_4, colData = Non_GSE234729_Samples, design = ~ disease_group + predicted_fetal_sex)

# pre_filter <- rowSums(counts(RNA_WG_DESeq_4)) >=10
# RNA_WG_DESeq_4 <- RNA_WG_DESeq_4[pre_filter,]
# RNA_Lip_DESeq_4 <- RNA_WG_DESeq_4[rownames(RNA_WG_DESeq_4) %in% lipid_ensembl_list$gene_id]

# RNA_Lip_DESeq_4$disease_group <- relevel(RNA_Lip_DESeq_4$disease_group, ref = "Control")
# RNA_Lip_DESeq_4 <- DESeq(RNA_Lip_DESeq_4)

# vsd_Lip_RNA_batch_uncorrected_4 <- vst(RNA_Lip_DESeq_4, blind = FALSE)

# PCA_batch_uncorrected_4_study <- plotPCA(vsd_Lip_RNA_batch_uncorrected_4, intgroup="GSE_number")

# png("./PCA_batch_uncorrected_4_study.png", height = 9, width = 15, units = "in", res = 300)
# grid.arrange(PCA_batch_uncorrected_4_study, nrow = 1)
# dev.off()

# PCA_batch_uncorrected_4_disease_group <- plotPCA(vsd_Lip_RNA_batch_uncorrected_4, intgroup = "disease_group")

# png("./PCA_batch_uncorrected_all4_disease_group.png", height = 9, width = 15, units = "in", res = 300)
# grid.arrange(PCA_batch_uncorrected_4_disease_group, nrow = 1)
# dev.off()

# #After Controlling for Batch Effects 
# vst_Lip_RNA_batch_uncorrected_4_matrix <- assay(vsd_Lip_RNA_batch_uncorrected_4)

# #vst_Lip_RNA_batch_corrected_4 <- limma::removeBatchEffect(vst_Lip_RNA_batch_uncorrected_4_matrix, batch=vsd_Lip_RNA_batch_uncorrected_4$GSE_number, batch2=vsd_Lip_RNA_batch_uncorrected_4$predicted_fetal_sex, design = model.matrix(~ disease_group, colData(vsd_Lip_RNA_batch_uncorrected_4)))
# vst_Lip_RNA_batch_corrected_4 <- limma::removeBatchEffect(vst_Lip_RNA_batch_uncorrected_4_matrix, batch=vsd_Lip_RNA_batch_uncorrected_4$GSE_number, batch2=vsd_Lip_RNA_batch_uncorrected_4$predicted_fetal_sex, design = model.matrix(~ vsd_Lip_RNA_batch_uncorrected_4$disease_group))


# vsd_Lip_RNA_batch_corrected_4 <- vsd_Lip_RNA_batch_uncorrected_4

# assay(vsd_Lip_RNA_batch_corrected_4) <- vst_Lip_RNA_batch_corrected_4

# PCA_batch_corrected_4_study <- plotPCA(vsd_Lip_RNA_batch_corrected_4, intgroup="GSE_number")

# png("./PCA_batch_corrected_4_study.png", height = 9, width = 15, units = "in", res = 300)
# grid.arrange(PCA_batch_corrected_4_study, nrow = 1)
# dev.off()

# PCA_batch_corrected_4_fetal_sex <- plotPCA(vsd_Lip_RNA_batch_corrected_4, intgroup = "predicted_fetal_sex")

# png("./PCA_batch_corrected_4_fetalsex.png", height = 9, width = 15, units = "in", res = 300)
# grid.arrange(PCA_batch_corrected_4_fetal_sex, nrow = 1)
# dev.off()

# PCA_batch_corrected_4_disease_group <- plotPCA(vsd_Lip_RNA_batch_corrected_4, intgroup = "disease_group")

# png("./PCA_batch_corrected_4_disease_group.png", height = 9, width = 15, units = "in", res = 300)
# grid.arrange(PCA_batch_corrected_4_disease_group, nrow = 1)
# dev.off()

### GSE218039
gse218039 <- getGEO('GSE218039', GSEMatrix = TRUE)
titles_gse218039 <- pData(phenoData(gse218039[[1]]))[1:3,] 
metadata218039 <- pData(phenoData(gse218039[[1]]))[,c("title","geo_accession","source_name_ch1", "molecule_ch1", "platform_id", "instrument_model", "library_strategy", "cell type:ch1", "preeclampsia:ch1", "rop:ch1", "Sex:ch1")]
metadata218039$GSE_number <- "GSE218039"
write.csv(metadata218039, "./metadata218039.csv")








# download raw data from: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE75010&format=file
# transfer to server

# Code adapted from "An end to end workflow for differential gene expression using Affymetrix microarrays", (Klaus & Reisenauer 2023)

BiocManager::install("maEndToEnd")

#General Bioconductor packages
    library(Biobase)
    library(oligoClasses)
     
#Annotation and data import packages
    library(ArrayExpress)
    library(pd.hugene.1.0.st.v1)
    library(hugene10sttranscriptcluster.db)
     
#Quality control and pre-processing packages
    library(oligo)
    library(arrayQualityMetrics)
    library(gridExtra) #version 2.3
    library(ggrepel) 

     
#Analysis and statistics packages
    library(limma)
    library(topGO)
    library(ReactomePA)
    library(clusterProfiler)
     
#Plotting and color options packages
    library(gplots)
    library(ggplot2)
    library(geneplotter)
    library(RColorBrewer)
    library(pheatmap)
    library(enrichplot)
     
#Formatting/documentation packages
   #library(rmarkdown)
   #library(BiocStyle)
    library(dplyr)
    library(tidyr)

#Helpers:
    library(stringr)
    library(matrixStats)
    library(genefilter)
    library(openxlsx)
   #library(devtools)


#Step 1: Download Dataset
library(GEOquery)
gse75010 <- getGEO('GSE75010', GSEMatrix = TRUE)

# identify the columns you want to extract for your metadata sheet
titles_gse75010 <- pData(phenoData(gse75010[[1]]))[1:3,]
titles_gse75010

#this is the Sample and Data Relationship Format (SDRF), phenoData
ga_data <- pData(phenoData(gse75010[[1]]))[,c("title", "ga (day):ch1","ga (week):ch1")]
ga_data$seven <- 7
names(ga_data)[names(ga_data) == "ga (day):ch1"] <- "ga_day"
names(ga_data)[names(ga_data) == "ga (week):ch1"] <- "ga_week"
ga_data$ga_day <- as.numeric(ga_data$ga_day)
ga_data$ga_week <- as.numeric(ga_data$ga_week)
ga_data$day_decimal<- ga_data$ga_day / 7
ga_data$gestational_age <- ga_data$ga_week + ga_data$day_decimal

metadata75010 <- pData(phenoData(gse75010[[1]]))[,c("title","geo_accession", "apgar score (1 min):ch1", "apgar score (5 min):ch1", "attempted vaginal delivery:ch1", "chorioamnionitis diagnosis:ch1", "molecule_ch1", "platform_id", "diagnosis:ch1", "hellp diagnosis:ch1", "infant gender:ch1","iugr diagnosis:ch1", "maternal age:ch1", "maternal blood type:ch1", "maternal bmi:ch1", "maternal ethnicity:ch1", "maximum diastolic bp:ch1", "maximum systolic bp:ch1", "mean umbilical pi:ch1", "mean uterine pi:ch1", "mode of delivery:ch1", "mode proteinuria:ch1", "newborn weight z-score:ch1", "nicu transfer:ch1", "placental weight z-score:ch1", "previous hypertensive pregnancy:ch1","previous miscarriage:ch1", "previous nulliparity:ch1", "umbilical cord diameter:ch1", "tissue:ch1", "source_name_ch1")]
metadata75010$gestational_age <- ga_data$gestational_age 
write.csv(metadata75010, "./metadata75010.csv")

#manually added column called "file_name" in metadata sheet with each sample's corresponding .CEL file name 
SDRF <- read.csv("./metadata75010.csv")
row.names(SDRF) <- SDRF$X

raw_data <- oligo::read.celfiles(filenames=file.path("./raw_data", SDRF$file_path), verbose = TRUE)
phenoData <- SDRF
row.names(phenoData) <- sampleNames(raw_data)
pd <- AnnotatedDataFrame(data = phenoData)
phenoData(raw_data) <- pd

#Quality Control

Biobase::exprs(raw_data)[1:5, 1:5]
#rows represent microarray probes, columns represent one microarrary (one sample)

exp_raw <- log2(Biobase::exprs(raw_data))
PCA_raw <- prcomp(t(exp_raw))

percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2],
                    Disease = pData(raw_data)$diagnosis.ch1,
                    Individual = pData(raw_data)$geo_accession,
                    Sex = pData(raw_data)$infant.gender.ch1,
                    BMI = pData(raw_data)$maternal.bmi.ch1,
                    Ethnicity = pData(raw_data)$maternal.ethnicity.ch1
)


PCA_disease <- ggplot(dataGG, aes(PC1, PC2)) +
geom_point(aes(color = Disease)) +
  ggtitle("PCA plot of the log-transformed raw expression data") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_fill_manual(values = c("#8a00c4","#d02670"))


png("./PCA_Affymetrix_Lip_PE_disease.png", height = 10, width = 10, units = "in", res = 300)
grid.arrange(PCA_disease, nrow = 1)
dev.off()

ggsave(PCA_disease, file = "./PCA_Affymetrix_Lip_PE_disease.png")

#Boxplot of log2-intensities 

png("./Boxplot_Affymetrix_Lip_PE_disease.png", height = 10, width = 10, units = "in", res = 300)
oligo::boxplot(exp_raw, target = "core", 
               main = "Boxplot of log2-intensitites for the raw data"
)
dev.off()

## Quality Control Metrics

arrayQuality <- arrayQualityMetrics(expressionset = raw_data,
    do.logtransform = TRUE,
    intgroup = c("diagnosis.ch1")
)

## Based on arrayQualityMetrics, need to remove array numbers 3, 35, 36, 56, 65, 68, 78, 85, 116, 141, 142 

outlier_samples <- c("GSM1940494_10119_HuGene-1_0-st-v1_.CEL", "GSM1940526_13700_HuGene-1_0-st-v1_.CEL", "GSM1940527_13746_HuGene-1_0-st-v1_.CEL", "GSM1940547_15015_HuGene-1_0-st-v1_.CEL", "GSM1940556_15811_HuGene-1_0-st-v1_.CEL", "GSM1940559_15873_HuGene-1_0-st-v1_.CEL", "GSM1940569_4520_HuGene-1_0-st-v1_.CEL", "GSM1940576_5157_HuGene-1_0-st-v1_.CEL", "GSM1940607_7668_HuGene-1_0-st-v1_.CEL", "GSM1940632_8802_HuGene-1_0-st-v1_.CEL", "GSM1940633_8971_HuGene-1_0-st-v1_.CEL")
filtered_SDRF <- SDRF[!SDRF$file_path %in% outlier_samples,]

filtered_raw_data <- oligo::read.celfiles(filenames=file.path("./raw_data", filtered_SDRF$file_path), verbose = TRUE)
filtered_phenoData <- filtered_SDRF
row.names(filtered_phenoData) <- sampleNames(filtered_raw_data)
filtered_pd <- AnnotatedDataFrame(data = filtered_phenoData)
phenoData(filtered_raw_data) <- filtered_pd

# Relative log expression data quality analysis
affymetrix_set <- oligo::rma(filtered_raw_data, target = "core", normalize = FALSE)

row_medians_assayData <- 
  Biobase::rowMedians(as.matrix(Biobase::exprs(affymetrix_set))
)

RLE_data <- sweep(Biobase::exprs(affymetrix_set), 1, row_medians_assayData)

RLE_data <- as.data.frame(RLE_data)
RLE_data_gathered <- 
  tidyr::gather(RLE_data, patient_array, log2_expression_deviation)

png("./RLE_Affymetrix_Lip_PE.png", height = 10, width = 40, units = "in", res = 300)
ggplot2::ggplot(RLE_data_gathered, aes(patient_array,
                                       log2_expression_deviation)) + 
  geom_boxplot(outlier.shape = NA) + 
  ylim(c(-2, 2)) + 
  theme(axis.text.x = element_text(colour = "purple", 
                                  angle = 60, size = 6.5, hjust = 1 ,
                                  face = "bold"))
dev.off()

#RMA calibration (background-correction, normalization, summarization - translates probe measurements to genes)

affymetrix_set_norm <- oligo::rma(filtered_raw_data, target = "core")

exp_affymetrix_Lip_PE <- Biobase::exprs(affymetrix_set_norm)
PCA <- prcomp(t(exp_affymetrix_Lip_PE), scale = FALSE)
percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
                    Disease = 
                     Biobase::pData(affymetrix_set_norm)$diagnosis.ch1
)

png("./PCA_Affymetrix_Lip_PE_Norm.png", height = 10, width = 10, units = "in", res = 300)

ggplot(dataGG, aes(PC1, PC2)) +
      geom_point(aes(colour = Disease)) +
  ggtitle("PCA plot of the calibrated, summarized data") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = sd_ratio) +
  scale_color_manual(values = c("darkorange2", "dodgerblue4"))
dev.off()


#Filtering based on intensity (filtering lowly expressed genes - low median expression)

affymetrix_Lip_PE_medians <- rowMedians(Biobase::exprs(affymetrix_set_norm))

##visually assess where the first peak on the histogram ends (that represents the background intensities)
png("./Histogram_Med_Intensities_Lip_PE.png", height = 10, width = 10, units = "in", res = 300)
hist_res <- hist(affymetrix_Lip_PE_medians, 100, col = "cornsilk1", freq = FALSE, main = "Histogram of the median intensities", border = "antiquewhite4", xlab = "Median intensities")
hist_res
dev.off()

man_threshold <- 4
hist_res <- hist(affymetrix_Lip_PE_medians, 100, col = "cornsilk1", freq = FALSE, main = "Histogram of the median intensities", border = "antiquewhite4", xlab = "Median intensities")
abline(v = man_threshold, col = "coral4", lwd = 2)

no_of_samples <- table(paste0(pData(affymetrix_set_norm)$diagnosis.ch1))
##non-PE = 72, PE = 74

##requires gene expression levels to be above noise threshold in at least 72 samples --> ensures that if a gene is only expressed in that minimum number of samples, it is retained
samples_cutoff <- min(no_of_samples)
idx_man_threshold <- apply(Biobase::exprs(affymetrix_set_norm), 1, function(x){sum(x > man_threshold) >= samples_cutoff})
table(idx_man_threshold) #FALSE: 2589 genes, TRUE: 30708 genes
affymetrix_Lip_PE_manfiltered <- subset(affymetrix_set_norm, idx_man_threshold)

#Annotation of transcript clusters
anno_affymetrix <- AnnotationDbi::select(hugene10sttranscriptcluster.db,
                                  keys = (featureNames(affymetrix_Lip_PE_manfiltered)),
                                  columns = c("SYMBOL", "GENENAME","CHR", "ENSEMBL", "PMID"),
                                  keytype = "PROBEID"
) #dim 3552715 6 
anno_affymetrix <- subset(anno_affymetrix, !is.na(ENSEMBL)) #dim 3534693 6
anno_grouped <- group_by(anno_affymetrix, PROBEID)
anno_summarized <- 
  dplyr::summarize(anno_grouped, no_of_matches = n_distinct(ENSEMBL)
)
dim(anno_summarized) #dim: 21271 2
anno_filtered <- filter(anno_summarized, no_of_matches > 1)
dim(anno_filtered) #dim: 3335 2
probe_stats <- anno_filtered 
final_probes <- anno_summarized %>% dplyr::filter(!PROBEID %in% probe_stats$PROBEID)
dim(final_probes) #dim: 17936 2

nrow(probe_stats) #3335
nrow(final_probes) #17936
ids_to_include <- (featureNames(affymetrix_Lip_PE_manfiltered) %in% final_probes$PROBEID)

table(ids_to_include) #FALSE: 12772 ENSEMBL IDs, TRUE: 17936 ENSEMBL IDs

affymetrix_Lip_PE_final <- subset(affymetrix_Lip_PE_manfiltered, ids_to_include)
validObject(affymetrix_Lip_PE_final) #TRUE

##exclude those same genes from the annotation file
head(anno_affymetrix)
anno_affymetrix_distinct <- anno_affymetrix[anno_affymetrix$PROBEID %in% anno_summarized$PROBEID,]
anno_affymetrix_distinct <- anno_affymetrix_distinct %>%
    dplyr::distinct(PROBEID, .keep_all = TRUE
)
fData(affymetrix_Lip_PE_final)$PROBEID <- rownames(fData(affymetrix_Lip_PE_final))
fData(affymetrix_Lip_PE_final) <- left_join(fData(affymetrix_Lip_PE_final), anno_affymetrix_distinct)
rownames(fData(affymetrix_Lip_PE_final)) <- fData(affymetrix_Lip_PE_final)$PROBEID 
validObject(affymetrix_Lip_PE_final) #dim 17936 6







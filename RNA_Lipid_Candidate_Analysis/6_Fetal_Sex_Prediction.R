#Code adapted from Tanya Phung's github code: SexChrLab/SexInference 
#Importing read counts from feature_counts output

BiocManager::install("tximport")
BiocManager::install("GenomicFeatures")
BiocManager::install("tidyr")

library(tximport)
library(GenomicFeatures)
library(tidyr)

gtf_file = "/workspace/lab/wilsonslab/datalake-wilsonslab/2025_RNA_Lipid_Candidate/genome_mapping/gencode.v48.primary_assembly.annotation.gtf"
samples_name = "./sample_name_SexInference.csv"
counts_file = "../read_counts.txt"
outfile = "./data_for_regression_2025.csv"

txdb = makeTxDbFromGFF(gtf_file)

k = keys(txdb, keytype = "TXNAME")

# write output with gene ID and transcript ID   
tx2gene = select(txdb, k, "GENEID", "TXNAME") 

samples_df = read.csv(samples_name)

featureCounts_data <- read.table(counts_file, header = TRUE, row.names = 1, sep = "\t")
counts_matrix <- as.matrix(featureCounts_data[3:nrow(featureCounts_data), c(6:ncol(featureCounts_data))])

counts_df = as.data.frame(counts_matrix)
counts_df$gene_id = rownames(counts_df)

XIST = subset(counts_df, counts_df$gene_id=="ENSG00000229807.15")
EIF1AY = subset(counts_df, counts_df$gene_id=="ENSG00000198692.10")
KDM5D = subset(counts_df, counts_df$gene_id=="ENSG00000012817.16")
UTY = subset(counts_df, counts_df$gene_id=="ENSG00000183878.16")
DDX3Y = subset(counts_df, counts_df$gene_id=="ENSG00000067048.17")
RPS4Y1 = subset(counts_df, counts_df$gene_id=="ENSG00000129824.16")

nsample = nrow(samples_df)

XIST_long <- gather(XIST[1:nsample], factor_key=TRUE)
colnames(XIST_long) = c("sample_ids", "XIST")

EIF1AY_long <- gather(EIF1AY[1:nsample], factor_key=TRUE)
colnames(EIF1AY_long) = c("sample_ids", "EIF1AY")

KDM5D_long <- gather(KDM5D[1:nsample], factor_key=TRUE)
colnames(KDM5D_long) = c("sample_ids", "KDM5D")

UTY_long <- gather(UTY[1:nsample], factor_key=TRUE)
colnames(UTY_long) = c("sample_ids", "UTY")

DDX3Y_long <- gather(DDX3Y[1:nsample], factor_key=TRUE)
colnames(DDX3Y_long) = c("sample_ids", "DDX3Y")


RPS4Y1_long <- gather(RPS4Y1[1:nsample], factor_key=TRUE)
colnames(RPS4Y1_long) = c("sample_ids", "RPS4Y1")

sample_name <- read.csv("./sample_name_SexInference.csv")
almost_srr <- sub("../genome_mapping/markeddup_BAMs/","", sample_name$sample_name)
srr <- sub("_markdup.bam", "", almost_srr)

data = data.frame(srr, XIST_long$XIST, EIF1AY_long$EIF1AY, KDM5D_long$KDM5D, UTY_long$UTY, DDX3Y_long$DDX3Y, RPS4Y1_long$RPS4Y1, samples_df$sex)
colnames(data) = c("sample_name","XIST","EIF1AY","KDM5D","UTY","DDX3Y","RPS4Y1","sex")

write.table(data, outfile, quote = F, row.names = F, sep=",")

#----------------------------------------------------------------------------------------------------------------------------------------

#Sex_Inference_Model

ibrary(tidyverse)
library(dplyr)
library(caret)
library(glmnet)

training_model_data = "/workspace/lab/wilsonslab/datalake-wilsonslab/2025_RNA_Lipid_Candidate/read_counts/sex_prediction/SexInference/RNAseq/training_data/for_rna_sex_check.tsv"
experiment_data = "/workspace/lab/wilsonslab/datalake-wilsonslab/2025_RNA_Lipid_Candidate/read_counts/sex_prediction/data_for_regression_2025.csv"

# -----------------------------------
# Build the model using the GTEx data
# -----------------------------------
# Load the data and remove NAs
data = read.csv(training_model_data, sep="\t")

# Split the data into training and test set
set.seed(123)
training.samples <- data$sex %>% 
  createDataPartition(p = 0.8, list = FALSE)
train.data  <- data[training.samples, ]
test.data <- data[-training.samples, ]

# Dumy code categorical predictor variables
x <- model.matrix(sex~., train.data)[,-1]
# Convert the outcome (class) to a numerical variable
y <- ifelse(train.data$sex == "female", 1, 0)

cv.lasso <- cv.glmnet(x, y, alpha = 1, family = "binomial")
plot(cv.lasso)
cv.lasso$lambda.min

coef(cv.lasso, cv.lasso$lambda.min)
coef(cv.lasso, cv.lasso$lambda.1se)

# Final model with lambda.min
lasso.model <- glmnet(x, y, alpha = 1, family = "binomial",
                      lambda = cv.lasso$lambda.min)
# Make prediction on test data
x.test <- model.matrix(sex ~., test.data)[,-1]
probabilities <- lasso.model %>% predict(newx = x.test)
predicted.classes <- ifelse(probabilities > 0.5, "female", "male")
# Model accuracy
observed.classes <- test.data$sex
mean(predicted.classes == observed.classes)

# ----------------------
# Run on experiment data
# ----------------------
sample_chart <- read.csv("./sample_name_SexInference.csv")
almost_srr <- sub("../genome_mapping/markeddup_BAMs/","", sample_name$sample_name)
srr <- sub("_markdup.bam", "", almost_srr)

experiment_data <- read.csv(experiment_data)
sample_names <- experiment_data[, 1]
rownames(experiment_data) <- sample_names
experiment_data <- experiment_data[, -1]

# Make prediction on placenta data
x.experiment <- model.matrix(sex ~., experiment_data)[,-1]
probabilities <- lasso.model %>% predict(newx = x.experiment)
predicted.classes <- ifelse(probabilities > 0.5, "F", "M")

df_observed <- data.frame(sample_name = rownames(experiment_data), observed_sex = experiment_data$sex)
final_predicted_sex <- left_join(df_observed, df_predictions, by = "sample_name")
write.csv (final_predicted_sex, "predicted_sex_lip_RNA_Seq.csv")


# Model accuracy
observed.classes <- experiment_data$sex
mean(predicted.classes == observed.classes)


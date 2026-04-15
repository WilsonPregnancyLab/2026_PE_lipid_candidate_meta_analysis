library(tidyverse)
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





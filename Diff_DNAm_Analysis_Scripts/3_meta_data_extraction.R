# R Version - R4.4.1
# Packages
library(GEOquery) #version 2.72.0

# Extracting Metadata from GEO for each data set
# Metadata is qualitative information about each sample (ex. fetal sex, control/case, tissue type)

### GSE98224
gse98224 <- getGEO('GSE98224', GSEMatrix = TRUE)
## Preview metadata fields to chose which metadata you want to extract based on what is available
titles_gse98224 <- pData(phenoData(gse98224[[1]]))[1:3,] 
write.csv(titles_gse98224, "titles_gse98224.csv")
# Extract metadata
gse98224 <- getGEO('GSE98224', GSEMatrix = TRUE)
metadata98224 <- pData(phenoData(gse98224[[1]]))[,c("title","geo_accession","platform_id","bmi:ch1","diagnosis:ch1","ga day:ch1","ga week:ch1","gender:ch1","hellp:ch1","iugr:ch1", "maternal age:ch1","maternal ethnicity:ch1","max diastolic:ch1","max systolic:ch1","mod:ch1","mode proteinuria:ch1","previous hypertensive pregnancy:ch1","previous miscarriage:ch1","tissue:ch1")]
write.csv(metadata98224, "./metadata98224.csv")

### GSE100197
gse100197 <- getGEO('GSE100197', GSEMatrix = TRUE)
## Preview metadata fields
titles_gse100197 <- pData(phenoData(gse100197[[1]]))[1:3,] 
write.csv(titles_gse100197, "titles_gse100197.csv")
# Extract metadata
gse100197 <- getGEO('GSE100197', GSEMatrix = TRUE)
metadata100197 <- pData(phenoData(gse100197[[1]]))[,c("title","geo_accession","molecule_ch1","platform_id","fetal sex:ch1","gestational age:ch1","pathology group:ch1","sample tissue:ch1","subject id:ch1")]
write.csv(metadata100197, "./metadata100197.csv")

### GSE75196
gse75196 <- getGEO('GSE98224', GSEMatrix = TRUE)
## Preview metadata fields
titles_gse75196 <- pData(phenoData(gse75196[[1]]))[1:3,] 
write.csv(titles_gse75196, "titles_gse75196.csv")
#Extract metadata
gse75196 <- getGEO('GSE75196', GSEMatrix = TRUE)
metadata75196 <- pData(phenoData(gse75196[[1]]))[,c("title","geo_accession","platform_id","disease:ch1","gestation (wk):ch1","Sex:ch1","tissue:ch1")]
write.csv(metadata75196, "./metadata75196.csv")

### GSE125605
gse125605 <- getGEO('GSE125605', GSEMatrix = TRUE)
## Preview metadata fields
titles_gse125605 <- pData(phenoData(gse125605[[1]]))[1:3,] 
write.csv(titles_gse125605, "titles_gse125605.csv")
#Extract metadata
gse125605 <- getGEO('GSE125605', GSEMatrix = TRUE)
metadata125605 <- pData(phenoData(gse125605[[1]]))[,c("title","geo_accession","platform_id","source_name_ch1","description","gestational age (weeks+days):ch1")]
write.csv(metadata125605, "./metadata125605.csv")

# Illumina Sample Sheet
# Sample sheets are .csv files (1 per data set) that help link the metadata to each sample.
# This sample sheet will be the only non-idat file (a type of file with the methylation data) in your base directory (refer to fetal_sex_prediction.R)
# These sample sheets were made manually by copying over the metadata from the above files to a new .csv file.

## Each sample sheet per data set should have the following column names: 
# Sample_Name (ex. GSM3578100_7668610115_R01C02)
# Geo_Accession (ex. GSM3578100)
# Pathology_Group (ex. Control)
# Sentrix_ID (ex. 7668610115)
# Sentrix_Position (ex. R01C02)
# Fetal_Sex	(if available)
# Gestational_Age	(if available)
# Tissue (if available)
# *Other meta data columns you have available

## Manually inspect the sample sheet for NA values and make a new column titled 'Exclude' whose value is "Exclude" if they have NA in any of the columns, and "keep" if all information is present
# This will be used to exclude samples prior to our DNA methylation analysis
# Template is linked: 


# Metadata Sheet
# This is a .csv file that combines the metadata for all of the data sets into one file
# This file does not have to be saved in a specific location

## A metadata sheet should have the following information: 
# Sample_Name (ex. GSM3578100_7668610115_R01C02)
# GSE_number (ex. GSE75196)
# Geo_Accession (ex. GSM3578100)
# Platform (ex. GPL13534)
# Sentrix_ID (ex. 7668610115)
# Sentrix_Position (ex. R01C02)
# Pathology_Group (ex. Control)
# Fetal_Sex	(if available)
# Gestational_Age	(if available)
# Tissue (if available)
# *Other meta data columns you have available
# file saved as "Metadata_Sheet_lipid_preeclampsia_excluded_removed.csv"

# to check your work, the total number of rows should equal the total number of samples


















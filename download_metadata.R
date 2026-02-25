# R Version - R4.4.1
# Install Packages
BiocManager::install("GEOquery", update = F)
library(GEOquery) #version 2.72.0
library(dplyr) #version 1.1.4

# Extracting Metadata from GEO for each data set
# Metadata is qualitative information about each sample (ex. fetal sex, control/case, tissue type)
# Metadata Sheet is a .csv file that combines the metadata for all of the data sets into one file

### GSE279757
gse279757 <- getGEO('GSE279757', GSEMatrix = TRUE)
## Preview metadata fields to chose which metadata you want to extract based on what is available
titles_gse279757 <- pData(phenoData(gse279757[[1]]))[1:3,] 
# Create table with the chosen metadata fields for all the samples in the study
metadata279757 <- pData(phenoData(gse279757[[1]]))[,c("title","geo_accession","source_name_ch1","molecule_ch1","description","platform_id","instrument_model","library_strategy","tissue:ch1")]
# Add Column with GSE number 
metadata279757$GSE_number <- "GSE279757"
write.csv(metadata279757, "./metadata279757.csv")

#Close R
q()

#Go back into terminal and transfer metadata file onto your computer
scp eyerk@fhssuperdome.csu.mcmaster.ca:/workspace/lab/wilsonslab/eyerk \Users\kriem\Downloads\

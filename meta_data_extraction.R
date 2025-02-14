# R Version - R4.4.1
# Packages
library(GEOquery) 

# Extracting Metadata from GEO

# GSE98224
gse98224 <- getGEO('GSE98224', GSEMatrix = TRUE)
## Preview metadata fields
titles_gse98224 <- pData(phenoData(gse98224[[1]]))[1:3,] 
write.csv(titles_gse98224, "titles_gse98224.csv")
# Extract metadata
gse98224 <- getGEO('GSE98224', GSEMatrix = TRUE)
metadata98224 <- pData(phenoData(gse98224[[1]]))[,c("title","geo_accession","platform_id","bmi:ch1","diagnosis:ch1","ga day:ch1","ga week:ch1","gender:ch1","hellp:ch1","iugr:ch1", "maternal age:ch1","maternal ethnicity:ch1","max diastolic:ch1","max systolic:ch1","mod:ch1","mode proteinuria:ch1","previous hypertensive pregnancy:ch1","previous miscarriage:ch1","tissue:ch1")]
write.csv(metadata98224, "./metadata98224.csv")

# GSE100197
gse100197 <- getGEO('GSE100197', GSEMatrix = TRUE)
## Preview metadata fields
titles_gse100197 <- pData(phenoData(gse100197[[1]]))[1:3,] 
write.csv(titles_gse100197, "titles_gse100197.csv")
# Extract metadata
gse100197 <- getGEO('GSE100197', GSEMatrix = TRUE)
metadata100197 <- pData(phenoData(gse100197[[1]]))[,c("title","geo_accession","molecule_ch1","platform_id","fetal sex:ch1","gestational age:ch1","pathology group:ch1","sample tissue:ch1","subject id:ch1")]
write.csv(metadata100197, "./metadata100197.csv")

# GSE75196
gse75196 <- getGEO('GSE98224', GSEMatrix = TRUE)
## Preview metadata fields
titles_gse75196 <- pData(phenoData(gse75196[[1]]))[1:3,] 
write.csv(titles_gse75196, "titles_gse75196.csv")
#Extract metadata
gse75196 <- getGEO('GSE75196', GSEMatrix = TRUE)
metadata75196 <- pData(phenoData(gse75196[[1]]))[,c("title","geo_accession","platform_id","disease:ch1","gestation (wk):ch1","Sex:ch1","tissue:ch1")]
write.csv(metadata75196, "./metadata75196.csv")

# GSE125605
gse125605 <- getGEO('GSE125605', GSEMatrix = TRUE)
## Preview metadata fields
titles_gse125605 <- pData(phenoData(gse125605[[1]]))[1:3,] 
write.csv(titles_gse125605, "titles_gse125605.csv")
#Extract metadata
gse125605 <- getGEO('GSE125605', GSEMatrix = TRUE)
metadata125605 <- pData(phenoData(gse125605[[1]]))[,c("title","geo_accession","platform_id","source_name_ch1","description","gestational age (weeks+days):ch1")]
write.csv(metadata125605, "./metadata125605.csv")

# Create 2-way frequency tables

# GSE98224
GSE_98224 <- read.csv("metadata98224.csv")
data_98224 <- table (GSE_98224$diagnosis, GSE_98224$gender)
write.csv (data_98224, "2way_98224.csv")

# GSE75196
GSE_75196 <- read.csv("metadata75196.csv")
data_75196 <- table (GSE_75196$disease, GSE_75196$Sex)
write.csv (data_75196, "2way_75196.csv")

# GSE100197
GSE_100197 <- read.csv("metadata100197.csv")
data_100197 <- table (GSE_100197$pathology_group, GSE_100197$fetal_sex)
write.csv (data_100197, "2way_100197.csv")

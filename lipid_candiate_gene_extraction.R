# downloaded datasets from https://maayanlab.cloud/Harmonizome/dataset/ --> GO Biological Process Annotations 2023, GO Cellular Component Annotations 2023  GAD Gene-Disease Associations, Reactome Pathways
# .txt file names were "gene_attribute_edges_GO2023_biological_processes.txt", "gene_attribute_edges_GO2023_cellular_components.txt", "gene_attribute_edges_GAD", "gene_attribute_edges_Reactome_Pathways"
# Based on literature review, saw some pathways I’m interested in looking at: 
# PPAR-γ, FFA, cholesterol, triglycerides, phospholipids, phosphatidylcholine, lysophosphatidylcholine, sphingolipids, ceramides, lipoxygenases, lipid peroxides

gene_GAD <- read.delim("gene_attribute_edges_GAD.txt") 
lipid_gene_GAD <- gene_GAD[grep("triglyc|phospholipid|sphingo|phospha|cholesterol|lipid|ceramide|fatty acid|triacyl|perioxisome",gene_GAD$Disease, ignore.case = TRUE), ]
write.csv(lipid_gene_GAD,"lipid_gene_GAD.csv") 

# do the same for the 3 other .txt files
# after this, opened the csv. files manually looked at 'GeneID' and 'GeneSymbol' columns and changed those names to format as 'Gene.ID' and 'Gene.Symbol' so that they would match between files
# named as “lipid_gene_GAD_1.csv”, etc.

# define variables for file merging

GAD <- read.csv("lipid_gene_GAD_1.csv")
GO_BP <-read.csv("lipid_gene_GO_BP_1.csv")
GO_CC <-read.csv("lipid_gene_GO_CC_1.csv")
Reactome <- read.csv("lipid_gene_Reactome_1.csv")
lipids_combined_files <- list (GO_BP,GO_CC,Reactome, GAD)
lipid_merge <- Reduce(function(x, y) merge(x, y, all=TRUE), lipids_combined_files)
write.csv(lipid_merge,"lipid_merge.csv")

# remove duplicates based on 'Gene.ID'

lipid_merge <- read.csv(“lipid_merge.csv”)
lipid_genes_unique <- lipid_merge[!duplicated(lipid_merge$Gene.ID), ]
write.csv(lipid_genes_unique,"lipid_genes_unique.csv")

# created new file called “paper_lipid_genes.csv” and manually added unique Genes Symbols from lit review to the bottom of the “lipid_gene_unique.csv” file.

# isolate genes associated with preeclampsia comorbidities

gene_GAD <- read.delim("gene_attribute_edges_GAD.txt")
comorb_gene_GAD <- gene_GAD[grep("obesity|diabete|hypertension",gene_GAD$Disease, ignore.case = TRUE), ]
write.csv (comorb_gene_GAD,"comorb_gene_GAD.csv")
comorb_gene_GAD <- read.csv("comorb_gene_GAD.csv")
comorb_genes_unique <- comorb_gene_GAD[!duplicated(comorb_gene_GAD$Gene.ID), ]
write.csv (comorb_genes_unique,"comorb_genes_unique.csv")

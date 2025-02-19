# downloaded datasets from https://maayanlab.cloud/Harmonizome/dataset/ --> GO Biological Process Annotations 2023 (https://maayanlab.cloud/Harmonizome/dataset/GO+Biological+Process+Annotations+2023), GO Cellular Component Annotations 2023 (https://maayanlab.cloud/Harmonizome/dataset/GO+Cellular+Component+Annotations+2023), 
  GAD Gene-Disease Associations (https://maayanlab.cloud/Harmonizome/dataset/GAD+Gene-Disease+Associations), Reactome Pathways 2024
# 'Gene-Attribute Edge list’ file for each data set was downloaded.
# .txt file names were "gene_attribute_edges_GO2023_biological_processes.txt", "gene_attribute_edges_GO2023_cellular_components.txt", "gene_attribute_edges_GAD", "gene_attribute_edges_Reactome_Pathways"
# Based on literature review looking lipid abnormalities in preeclampsia patients, saw some pathways I’m interested in looking at: 
# PPAR-γ, free fatty acid, cholesterol (LDL & HDL), triglycerides, phospholipids, phosphatidylcholine, lysophosphatidylcholine, sphingolipids, ceramide, lipoxygenases, lipid peroxides, peroxisome, triacylglycerols. 

# define the txt file as a variable so that we can manipulate it
gene_GO_BP <- read.delim("gene_attribute_edges_GO2023_biological_processes.txt") 
# use the grep function to identify any of these terms separated by "|" (meaning "or") in the GO.Biological.Process column of gene_GO_BP. We don't care for the case.
lipid_gene_GO_BP <- gene_GO_BP[grep("triglyc|phospholipid|sphingo|phospha|cholesterol|lipid|ceramide|fatty acid|triacyl|perioxisome",gene_GO_BP$GO.Biological.Process, ignore.case = TRUE), ]
# write a new csv file that has only our genes whose descriptions are filtered for our terms 
write.csv(lipid_gene_GO_BP,"lipid_gene_GO_BP.csv") 

# repeat the above steps for each of the 4 data sets
gene_GO_CC <- read.delim("gene_attribute_edges_GO2023_cellular_components.txt") 
lipid_gene_GO_CC <- gene_GO_CC[grep("triglyc|phospholipid|sphingo|phospha|cholesterol|lipid|ceramide|fatty acid|triacyl|perioxisome",gene_GO_CC$GO.Cellular.Component, ignore.case = TRUE), ]
write.csv(lipid_gene_GO_CC,"lipid_gene_GO_CC.csv") 

# here, we have 1 extra step where we rename the GeneSym and GeneID columns to match those from the previous 2 datas sets. this allows for an easier time merging the data sets later on. 
gene_GAD <- read.delim("gene_attribute_edges_GAD.txt") 
colnames(gene_GAD)[colnames(gene_GAD) == "GeneSym"] <- "Gene.Symbol"
colnames(gene_GAD)[colnames(gene_GAD) == "GeneID"] <- "Gene.ID"
lipid_gene_GAD <- gene_GAD[grep("triglyc|phospholipid|sphingo|phospha|cholesterol|lipid|ceramide|fatty acid|triacyl|perioxisome",gene_GAD$Disease, ignore.case = TRUE), ]
write.csv(lipid_gene_GAD,"lipid_gene_GAD.csv") 

gene_Reactome <- read.delim("gene_attribute_edges_Reactome_Pathways.txt") 
colnames(gene_Reactome)[colnames(gene_Reactome) == "GeneSym"] <- "Gene.Symbol"
colnames(gene_Reactome)[colnames(gene_Reactome) == "GeneID"] <- "Gene.ID"
lipid_gene_Reactome <- gene_Reactome[grep("triglyc|phospholipid|sphingo|phospha|cholesterol|lipid|ceramide|fatty acid|triacyl|perioxisome",gene_Reactome$Pathway, ignore.case = TRUE), ]
write.csv(lipid_gene_Reactome,"lipid_gene_Reactome.csv") 

# define variables for file merging

GAD <- read.csv("lipid_gene_GAD.csv")
GO_BP <-read.csv("lipid_gene_GO_BP.csv")
GO_CC <-read.csv("lipid_gene_GO_CC.csv")
Reactome <- read.csv("lipid_gene_Reactome.csv")

# merge the lipid files
lipids_combined_files <- list (GO_BP,GO_CC,Reactome, GAD)
lipid_merge <- Reduce(function(x, y) merge(x, y, all=TRUE), lipids_combined_files)
write.csv(lipid_merge,"lipid_merge.csv")

# remove duplicates based on 'Gene.ID'. This will isolate the unique genes based on Gene.ID.

lipid_merge <- read.csv(“lipid_merge.csv”)
lipid_genes_unique <- lipid_merge[!duplicated(lipid_merge$Gene.ID), ]
write.csv(lipid_genes_unique,"lipid_genes_unique.csv")

# created new file called “paper_lipid_genes.csv” and manually added unique Genes Symbols from lit review to the bottom of the “lipid_gene_unique.csv” file. 
# references noted in thesis report 
                      
# isolate genes associated with preeclampsia comorbidities - similar to previous extraction however only use GAD data set

gene_GAD <- read.delim("gene_attribute_edges_GAD.txt")
comorb_gene_GAD <- gene_GAD[grep("obesity|diabete|hypertension",gene_GAD$Disease, ignore.case = TRUE), ]
write.csv (comorb_gene_GAD,"comorb_gene_GAD.csv")
comorb_gene_GAD <- read.csv("comorb_gene_GAD.csv")
comorb_genes_unique <- comorb_gene_GAD[!duplicated(comorb_gene_GAD$Gene.ID), ]
write.csv (comorb_genes_unique,"comorb_genes_unique.csv")

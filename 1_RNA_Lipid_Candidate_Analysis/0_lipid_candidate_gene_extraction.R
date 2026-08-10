##
# downloaded datasets from https://maayanlab.cloud/Harmonizome/dataset/ --> GAD Gene-Disease Associations (https://maayanlab.cloud/Harmonizome/dataset/GAD+Gene-Disease+Associations), Reactome Pathways 2024 (https://maayanlab.cloud/Harmonizome/dataset/Reactome+Pathways+2024)
# 'Gene-Attribute Edge list’ file for each data set was downloaded.
# .txt file names were "gene_attribute_edges_GAD", "gene_attribute_edges_Reactome_Pathways"
# Based on literature review looking lipid abnormalities in preeclampsia patients, saw some pathways I’m interested in looking at: 
# PPAR-γ, free fatty acid, cholesterol (LDL & HDL), triglycerides, phospholipids, phosphatidylcholine, lysophosphatidylcholine, sphingolipids, ceramide, lipoxygenases, lipid peroxides, peroxisome, triacylglycerols. 

# define the txt file as a variable so that we can manipulate it
gene_GAD <- read.delim("gene_attribute_edges_GAD.txt") 
# change name of "Gene Symbol" column to match format of other files
colnames(gene_GAD)[colnames(gene_GAD) == "GeneSym"] <- "gene_symbol"
# use the grep function to identify any of these terms separated by "|" (meaning "or") in the GO.Biological.Process column of gene_GO_BP. We don't care for the case.
lipid_gene_GAD <- gene_GAD[grep("triglyc|phospholipid|sphingo|cholesterol|lipid|lipoprotein|ceramide|fatty acid|triacylgly|peroxisome",gene_GAD$Disease, ignore.case = TRUE), ]
# write a new csv file that has only our genes whose descriptions are filtered for our terms 
write.csv(lipid_gene_GAD,"lipid_gene_GAD.csv") 

# repeat for Reactome database
gene_Reactome <- read.delim("gene_attribute_edges_Reactome_Pathways.txt") 
colnames(gene_Reactome)[colnames(gene_Reactome) == "GeneSym"] <- "gene_symbol"
lipid_gene_Reactome <- gene_Reactome[grep("triglyc|phospholipid|sphingo|cholesterol|lipid|lipoprotein|ceramide|fatty acid|triacylgly|peroxisome",gene_Reactome$Pathway, ignore.case = TRUE), ]
write.csv(lipid_gene_Reactome,"lipid_gene_Reactome.csv") 

# define variables for file merging

GAD <- read.csv("lipid_gene_GAD.csv")
Reactome <- read.csv("lipid_gene_Reactome.csv")
# downloaded lipid-related gene list from LIPID MAPS Gene/Proteome Database: https://lipidmaps.org/databases/lmpd/download
LMPD <- read.csv("LMPD_040215.csv")
human_LMPD <- subset(LMPD, LMPD$species == "Human")
# manually exported from lipid-related gene list DBLiPro Lipoprotein database: http://lipid.cloudna.cn/lipoprotein 
DBLiPro <- read.csv("DBLiPro_Lipoproteins.csv")
# searched through literature for lipid-related genes and gene lists 
Lit <- read.csv("Lipid_Related_Genes_Literature_Search.csv")

# merge the lipid files
lipid_gene_lists <- list (GAD, Reactome, human_LMPD, DBLiPro, Lit)
lipid_lists_merge <- Reduce(function(x, y) merge(x, y, all=TRUE), lipid_gene_lists)

# remove duplicates based on 'gene_symbol'. This will isolate the unique genes based on gene_symbol.

lipid_genes_unique <- lipid_lists_merge[!duplicated(lipid_lists_merge$gene_symbol), ]
condensed_lipid_genes_unique <- lipid_genes_unique[, (names(lipid_genes_unique) %in% c("gene_symbol", "Disease", "Pathway", "protein_name", "Protein.Names", "Title"))]
write.csv(condensed_lipid_genes_unique,"lipid_genes_unique.csv") #5510 unique lipid-related genes
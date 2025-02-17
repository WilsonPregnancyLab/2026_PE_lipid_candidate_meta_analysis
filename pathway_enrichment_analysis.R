# Part 1: NCBI Reference Subsetting
# the goal here is to add the NCBI annotations to the lipid genes (First, I'm going to do it for all of the lipid genes I tested - in whole_pop, M and F, then only for the significant ones in each group)

 if(!file.exists('Generic_human_ncbiIds_noParents.an.txt.gz')){
    system('wget https://gemma.msl.ubc.ca/annots/Generic_human_ncbiIds_noParents.an.txt.gz')}
NCBI <- read.table('Generic_human_ncbiIds_noParents.an.txt.gz', sep = '\t', header = T, quote = "")
placmet_wholepop_auto <- read.csv("placmet_wholepop_auto.csv")

all_wholepop_genes <- placmet_wholepop_auto 
all_wholepop_genes <- unique(all_wholepop_genes) #46613 lipid genes investigated

# Extracts the NCBI annotation file only for the genes that overlap our lipid genes. If the lipid gene is not present in the NCBI file, it will not be present in the dataframe. 
all_wholepop_reference <- NCBI[NCBI$GeneSymbols %in% all_wholepop_genes$Closest_TSS_gene_name, ] #2845 lipid genes w/ GO annots

write.csv(background_reference, "placmet_all_wholepop_NCBI_reference.csv")

# Part 2: Pathway Enrichment Analysis
devtools::install_github('PavlidisLab/ermineR')
library(ermineR)
library(dplyr)
library(ggplot2)
install.packages("rJava")

#I'M HAVING TROUBLE!!! JAVA ISN'T INSTALLING!!! I SHALL ASK KEATON
Sys.setenv('JAVA_HOME' = '/usr/lib/jvm/java-21-openjdk-amd64/bin/java')

all_wholepop_genes <- read.csv("placmet_wholepop_auto.csv")
all_wholepop_annotations <- read.csv("placmet_all_wholepop_NCBI_reference.csv")

all_wholepop_pr_out <- precRecall(annotation = all_wholepop_annotations, 
                    scores = all_wholepop_genes,
                    scoreColumn = 2,
                    logTrans = T,
                    bigIsBetter = F,
                    iterations = 10000)

all_wholepop_pathway_analysis <- as.data.frame(all_wholepop_pr_out$results)






















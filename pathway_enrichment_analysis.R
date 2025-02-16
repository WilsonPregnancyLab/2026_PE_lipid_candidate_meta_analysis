# Part 1: NCBI Reference Subsetting
# the goal here is to add the NCBI annotations to the lipid genes (here I'm doing it for wholepop labelled as bio sig)

 if(!file.exists('Generic_human_ncbiIds_noParents.an.txt.gz')){
    system('wget https://gemma.msl.ubc.ca/annots/Generic_human_ncbiIds_noParents.an.txt.gz')
NCBI <- read.table('Generic_human_ncbiIds_noParents.an.txt.gz', sep = '\t', header = T, quote = "")
placmet_wholepop_auto <- read.csv("placmet_wholepop_auto.csv")
background_genes <- placmet_wholepop_auto
background_reference <- NCBI[NCBI$GeneSymbols %in% background_genes$Closest_TSS_gene_name,]

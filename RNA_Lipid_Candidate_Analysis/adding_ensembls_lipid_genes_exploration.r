#Need to unearth the lipid related gene list to filter the read counts
lipid_gene_list <- read.csv("lipid_genes_unique_P.csv")

#Convert GENE_IDs to Ensembl IDs to identify gene identity on read_counts

library("AnnotationDbi")
library("org.Hs.eg.db")

#assess what options are available for keytypes
keytypes(org.Hs.eg.db)

#keys (what we want to lookup values for), columns (what we want returned back), keytype (the type of key we are passing)
lipid_gene_list$ENSEMBL <- mapIds(org.Hs.eg.db,
                    keys=lipid_gene_list$gene_symbol, 
                    column="ENSEMBL",
                    keytype="SYMBOL",
                    multiVals="list"
)


#biomaRt

BiocManager::install("biomaRt")
BiocManager::install("EnsDb.Hsapiens.v86")
BiocManager::install("rtracklayer")
BiocManager::install("data.table")

library(biomaRt)
library(EnsDb.Hsapiens.v86)
library(data.table)

lipid_gene_list <- read.csv("lipid_genes_unique_P.csv")
listEnsembl()
ensembl <- useEnsembl(biomart = "genes")
datasets <- listDatasets(ensembl)

ensembl_con <- useMart("ensembl", dataset = 'hsapiens_gene_ensembl')

LIPGENE_TO_ENSEMBL <- getBM(attributes = c('external_gene_name','ensembl_gene_id'),
        filters = "external_gene_name",
        values = lipid_gene_list$gene_symbol,
        mart = ensembl_con
)
write.csv(LIPGENE_TO_ENSEMBL, "LIPGENE_TO_ENSEMBL.csv")


read_counts <- read.csv("../read_counts/read_count_nohead.csv") #read_count file with no header (starts from column names)
filt_read_counts <- read_counts[2:nrow(read_counts), ]
filt_read_counts$Column1 <- sub("\\..*", "", filt_read_counts$Column1)


ENSEMBL_TO_ALLGENES <- getBM(attributes = c('ensembl_gene_id','external_gene_name'),
        filters = "ensembl_gene_id",
        values = filt_read_counts$Column1,
        mart = ensembl_con
)

write.csv(ENSEMBL_TO_ALLGENES, "ENSEMBL_TO_ALLGENES.csv")


#back to AnnotationDBI again

#Some symbols aren't matching to an ENSEMBL ID because they are aliases, these names will get me the aliases and hopefully extract the ENSEMBL IDs
official_to_alias <- mapIds(org.Hs.eg.db,
        keys = lipid_gene_list$gene_symbol,
        column = "SYMBOL",
        keytype = "ALIAS",
        multiVals = "list"
)
official_to_alias_list <- stack(official_to_alias)
write.csv(official_to_alias_list, "official_to_alias_list.csv")

#some weirdness with the annotation file, will try again with the gencode file i used to map the reads
library("rtracklayer")
library("dplyr")

gtf_file <- "../genome_mapping/gencode.v48.primary_assembly.annotation.gtf"
gtf_data <- import(gtf_file)
gtf_df <- as.data.frame(gtf_data)
gtf_df$gene_id <- sub("\\..*", "", gtf_df$gene_id)
gtf_df_dups_rem <- gtf_df[!duplicated(gtf_df$gene_id),]

read_counts <- read.csv("../read_counts/read_count_nohead.csv") #read_count file with no header (starts from column names)
colnames(read_counts) <- read_counts[1, ]
filt_read_counts <- read_counts[2:nrow(read_counts), ]
filt_read_counts$Geneid <- sub("\\..*", "", filt_read_counts$Geneid)
filt_read_counts <- merge(filt_read_counts,gtf_df_dups_rem[, c("gene_id", "gene_name")], join_by(Geneid == gene_id))

setDT(filt_read_counts)
setDT(gtf_df_dups_rem)
named_read_counts <- merge(x = filt_read_counts, y= gtf_df_dups_rem[, c("gene_id", "gene_name")], by.x = "Geneid", by.y = "gene_id")
ensembl_gene_names <- named_read_counts[,c("Geneid", "gene_name")]
write.csv(ensembl_gene_names, "ensembl_gene_names.csv")



ENSEMBL_TO_ALLGENES_gencode <- getBM(attributes = c('ensembl_gene_id','external_gene_name'),
        filters = "ensembl_gene_id",
        values = filt_read_counts$Column1,
        mart = ensembl_con
)

alias <- mapIds(org.Hs.eg.db,
        keys = lipid_gene_list$gene_symbol,
        column = "SYMBOL",
        keytype = "ALIAS",
        multiVals = "list"
)

alias_list <- stack(alias)

write.csv(alias_list, "alias_list.csv")

write.csv(ENSEMBL_TO_ALLGENES, "ENSEMBL_TO_ALLGENES.csv")



#mapping lipid_genes to gencode ensembl numbers

lipid_gene_list <- read.csv("lipid_genes_unique_P.csv")

setDT(lipid_gene_list)
setDT(gtf_df_dups_rem)
lipid_ensembl <- merge(x = lipid_gene_list, y= gtf_df_dups_rem[, c("gene_id", "gene_name")], by.x = "gene_symbol", by.y = "gene_name", all.x = TRUE)
#lipid_ensembl <- named_read_counts[,c("Geneid", "gene_name")]
write.csv(lipid_ensembl, "lipid_ensembl_gencode.csv")

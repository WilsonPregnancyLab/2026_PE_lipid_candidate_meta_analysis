library(rtracklayer)
library(data.table)
library(dplyr)
library(tidyverse)
library(DESeq2)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ggplot2) #version 3.5.1


gtf_file <- "./2025_RNA_Lipid_Candidate/genome_mapping/gencode.v48.primary_assembly.annotation.gtf"
gtf_data <- import(gtf_file)
gtf_df <- as.data.frame(gtf_data)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

auto_F_table <- read.csv("./Lip_PE_Met_RNA_auto_F.csv")
auto_M_table <- read.csv("./Lip_PE_Met_RNA_auto_M.csv")

WG_auto_F_table <- read.csv("./WG_PE_Met_RNA_auto_F.csv")
WG_auto_M_table <- read.csv("./WG_PE_Met_RNA_auto_M.csv")

placmet_M_fulldata_auto <- read.csv("/workspace/lab/wilsonslab/datalake-wilsonslab/PE_Lipid_Meta-analysis/2025_DNAm_Lipid_Candidate/Candidate_Diff_DNAm_Analysis_2025/placmet_M_fulldata_auto_rerun.csv")
wg_placmet_M_fulldata_auto <- read.csv("/workspace/lab/wilsonslab/datalake-wilsonslab/PE_Lipid_Meta-analysis/2025_DNAm_Lipid_Candidate/WG_Diff_DNAm_Analysis_2025/wg_placmet_M_fulldata_auto_rerun.csv")

# In F, no genes that are both significantly differentially expressed and methylated 
# Filter M tables to only include genes that are either trending towards or differentially expressed

include <- c("Increased_RNA_Expression", "Decreased_RNA_Expression", "Trending_Towards_Decreased_RNA_Expression", "Trending_Towards_Increased_RNA_Expression")
super_include <- c("Increased_RNA_Expression", "Decreased_RNA_Expression")

diff_exp_auto_M <- auto_M_table[auto_M_table$Expression_Status %in% super_include, ]
# diff_exp_WG_auto_M <- WG_auto_M_table[WG_auto_M_table$Expression_Status %in% super_include, ]


gene_gtf_df <- gtf_df[gtf_df$type == "gene", ]
gtf_df_cols <- c("seqnames", "start", "end", "type", "gene_id", "gene_name")
auto_M_gene <- gene_gtf_df[gene_gtf_df$gene_name %in% auto_M_table$gene_symbol, gtf_df_cols]


#CMIP test
auto_M_gene_CMIP <- auto_M_gene[auto_M_gene$gene_name == "CMIP", ]
CMIP_CpG <- placmet_M_fulldata_auto %>% filter(placmet_M_fulldata_auto$position %in% ((auto_M_gene_CMIP$start-3000):auto_M_gene_CMIP$end))
noDup_CMIP_CpG <- CMIP_CpG[!duplicated(CMIP_CpG$probes),]
sig_CMIP_CpG <- CMIP_CpG[CMIP_CpG$adj.P.Val < 0.05, ]

df_CMIP_all <- CMIP_CpG[,c("probes", "adj.P.Val", "deltaB", "region_overlap")]
df_CMIP_all$diffexp_gene <- "CMIP"

# df_CMIP_Dup <- CMIP_CpG[,c("probes", "deltaB")]
# df_CMIP <- df_CMIP_Dup[!duplicated(df_CMIP_Dup$probes),]
# df_CMIP <- as.data.frame(df_CMIP)
# row.names(df_CMIP) <- df_CMIP$probes
# df_CMIP$probes <- NULL
# heat_map_df_CMIP <- t(df_CMIP)

#heat_map_df_CMIP <- heat_map_df_CMIP[heat_map_df_CMIP$c("deltaB"),]

CMIP_data <- as.matrix(heat_map_df_CMIP)

png("./CMIP_heatmap.png")

ggplot(df_CMIP_all, aes(x = probes, y = diffexp_gene, fill = deltaB)) +
  geom_tile() +
  scale_fill_gradient2(low = "#8a00c4",
                       mid = "white",
                       high = "#d02670") +
    # coord_fixed() +
    guides(fill = guide_colourbar(barwidth = 0.5,
                                barheight = 20))
 


dev.off()



#function
deltaB_heatmap_full <- function(table, placmet_data) {
gene_gtf_df <- gtf_df[gtf_df$type == "gene", ]
gtf_df_cols <- c("seqnames", "start", "end", "type", "gene_id", "gene_name")
table <- table[!duplicated(table$gene_symbol),]
auto_gene <- gene_gtf_df[gene_gtf_df$gene_name %in% table$gene_symbol, gtf_df_cols]

CpG_list <- list()
for (i in 1:nrow(auto_gene)) {
    gene_start <- auto_gene[i, "start"] - 3000
    gene_end <- auto_gene[i, "end"]
    gene_name <- auto_gene[i, "gene_name"]

    CpGs_in_range <- placmet_data[placmet_data$position >= gene_start & placmet_data$position <= gene_end, ]
    CpGs_in_range <- cbind(CpGs_in_range, associated_gene = gene_name)

    CpG_list[[i]] <- CpGs_in_range

}

CpG_full_list <- do.call(rbind, CpG_list)

df_all_CpG <- CpG_full_list[,c("probes", "adj.P.Val", "deltaB", "region_overlap", "gene", "position", "associated_gene")]
df_all_CpG <- df_all_CpG[!duplicated(df_all_CpG$probes),]
df_all_CpG <- df_all_CpG %>% 
    group_by(associated_gene) %>% 
    mutate (number = row_number()) 
df_all_CpG <- as.data.frame(df_all_CpG)


return(df_all_CpG)

}

deltaB_hm_table <- deltaB_heatmap_full (auto_M_table, placmet_M_fulldata_auto)




deltaB_heatmap_around_sig <- function(table, placmet_data) {
gene_gtf_df <- gtf_df[gtf_df$type == "gene", ]
gtf_df_cols <- c("seqnames", "start", "end", "type", "gene_id", "gene_name")
table <- table[!duplicated(table$gene_symbol),]
auto_gene <- gene_gtf_df[gene_gtf_df$gene_name %in% table$gene_symbol, gtf_df_cols]
auto_gene <- merge(auto_gene, table[, c("gene_symbol","position")], by.x = "gene_name", by.y = "gene_symbol")

CpG_list <- list()
for (i in 1:nrow(auto_gene)) {
    gene_start <- auto_gene[i, "position"] - 1000
    gene_end <- auto_gene[i, "end"] + 1000
    gene_name <- auto_gene[i, "gene_name"]

    CpGs_in_range <- placmet_data[placmet_data$position >= gene_start & placmet_data$position <= gene_end, ]
    CpGs_in_range <- cbind(CpGs_in_range, associated_gene = gene_name)

    CpG_list[[i]] <- CpGs_in_range

}

CpG_full_list <- do.call(rbind, CpG_list)

df_all_CpG <- CpG_full_list[,c("probes", "adj.P.Val", "deltaB", "region_overlap", "gene", "position", "associated_gene")]
df_all_CpG <- df_all_CpG[!duplicated(df_all_CpG$probes),]

df_all_CpG <- df_all_CpG %>% 
    group_by(associated_gene) %>% 
    mutate (probe_number = row_number()) 
df_all_CpG <- as.data.frame(df_all_CpG)

return(df_all_CpG)
}

deltaB_hm_sig <- deltaB_heatmap_around_sig(auto_M_table, placmet_M_fulldata_auto)
GATA3 <- deltaB_hm_sig[deltaB_hm_sig$associated_gene == "GATA3", ]
CMIP <- deltaB_hm_sig[deltaB_hm_sig$associated_gene == "CMIP", ]
KIF26B <- deltaB_hm_sig[deltaB_hm_sig$associated_gene == "KIF26B", ]
C8orf58 <- deltaB_hm_sig[deltaB_hm_sig$associated_gene == "C8orf58", ]


png("./Lip_PE_CMIP_gene_region_sig.png")

ggplot(CMIP, aes(x = probes, y = associated_gene, fill = deltaB)) +
  geom_tile() +
  scale_fill_gradient2(low = "#8a00c4",
                       high = "#d02670") +
    #coord_fixed() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 2)) +
    guides(fill = guide_colourbar(barwidth = 0.5,
                                barheight = 5)) 
 


dev.off()


# Figuring Stuff Out
# gene_gtf_df <- gtf_df[gtf_df$type == "gene", ]
# gtf_df_cols <- c("seqnames", "start", "end", "type", "gene_id", "gene_name")
# auto_M_gene <- gene_gtf_df[gene_gtf_df$gene_name %in% auto_M_table$gene_symbol, gtf_df_cols]

# # auto_M_gene <- auto_M_gene[auto_M_gene$gene_name == "CMIP", ]

# CpG_list <- list()
# for (i in 1:nrow(auto_M_gene)) {
#     gene_start <- auto_M_gene[i, "start"] - 3000
#     gene_end <- auto_M_gene[i, "end"]
#     gene_name <- auto_M_gene[i, "gene_name"]

#     CpGs_in_range <- placmet_M_fulldata_auto[placmet_M_fulldata_auto$position >= gene_start & placmet_M_fulldata_auto$position <= gene_end, ]
#     CpGs_in_range <- cbind(CpGs_in_range, associated_gene = gene_name)

#     CpG_list[[i]] <- CpGs_in_range

# }

# CpG <- do.call(rbind, CpG_list)
# noDup_CpG <- CpG[!duplicated(CpG$probes),]
# sig_CpG <- CpG[CpG$adj.P.Val < 0.05, ]

# df_all_CpG <- CpG[,c("probes", "adj.P.Val", "deltaB", "region_overlap", "gene", "position", "associated_gene")]
# df_all_CpG <- df_all_CpG %>% 
#     group_by(associated_gene) %>% 
#     mutate (number = row_number()) 
# df_all_CpG <- as.data.frame(df_all_CpG)





png("./lipid_M_heatmap_number.png")

ggplot(df_all_CpG, aes(x = number, y = associated_gene, fill = deltaB)) +
  geom_tile() +
  scale_fill_gradient2(low = "#8a00c4",
                       high = "#d02670") +
    # coord_fixed() +
    guides(fill = guide_colourbar(barwidth = 0.5,
                                barheight = 20))
 


dev.off()
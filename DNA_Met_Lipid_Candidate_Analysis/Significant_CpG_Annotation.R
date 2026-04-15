#Adapted from mel-14's code Melanie_PhD_Code/IVF_MetaAnalysis/Significant_CpG_Annotation.R

BiocManager::install("IlluminaHumanMethylationEPICv2anno.20a1.hg38")
BiocManager::install("minfi")
BiocManager::install("annotatr")
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")

library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(minfi)
library(annotatr)
library(readr)
library(tibble)
library(tidyr)
library(tidyverse)

# R 4.5.1 - done locally 
# Load hg38 annotations for promoters and gene bodies
annotations <- c("hg38_genes_promoters", "hg38_genes_5UTRs","hg38_genes_3UTRs","hg38_genes_introns","hg38_basicgenes")
locations <- build_annotations(genome = "hg38", annotations = annotations)
locations_df <- as.data.frame(locations)


#R4.4.1 - done on Server 
cd /workspace/lab/wilsonslab/datalake-wilsonslab/2025_DNAm_Lipid_Candidate/Candidate_DNAm_Annotations

placmet_wholepop_auto_win <- read.csv("placmet_wholepop_auto_win.csv")
placmet_M_fulldata_auto   <- read.csv("placmet_M_fulldata_auto.csv")
placmet_F_fulldata_auto <- read.csv("placmet_F_fulldata_auto.csv")

supplementary2_Whole <- subset(placmet_wholepop_auto_win[placmet_wholepop_auto_win$adj.P.Val <0.05,])
supplementary2_Male   <- subset(placmet_M_fulldata_auto[placmet_M_fulldata_auto$adj.P.Val <0.05,])
supplementary2_Female <- subset(placmet_F_fulldata_auto[placmet_F_fulldata_auto$adj.P.Val <0.05,])

annotation19 <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
wholesig_19 <- as.data.frame(merge(supplementary2_Whole, annotation19, by.x = "probes", by.y = "Name", all.x = TRUE))
malesig_19 <- as.data.frame(merge(supplementary2_Male, annotation19, by.x = "probes", by.y = "Name", all.x = TRUE))
femalesig_19 <- as.data.frame(merge(supplementary2_Female, annotation19, by.x = "probes", by.y = "Name", all.x = TRUE))

#Function for making .bed files 
make_the_bed <- function(df = "df",
                         pos = "pos",
                         chr = "chr",
                         cpg = "cpg",
                         save_info = "/location/name.bed") {
  df2 <- df %>%
  mutate(posstart = pos-1) %>%
  mutate(posend = paste(
    pos
  )) %>%
  dplyr::select(chr, posstart, posend, cpg)

  write.table(df2,
    file = save_info,
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )
}

make_the_bed(wholesig_19, pos = "pos", chr = "chr", cpg = "probes", save_info = "./wholesig_19.bed")
make_the_bed(malesig_19, pos = "pos", chr = "chr", cpg = "probes", save_info = "./malesig_19.bed")
make_the_bed(femalesig_19, pos = "pos", chr = "chr", cpg = "probes", save_info = "./femalesig_19.bed")

unmake_the_bed <- function(file = "./name.bed" 
) {
  df <- read.table(file, sep = "\t")
  df %>%
    dplyr::select(
      chr = V1,
      pos = V3,
      CpG = V4
    )
}


#------------------------------------------------



whole_cpg_38pos <- unmake_the_bed(file = "./hglft_whole_cpg_38pos.bed")
female_cpg_38pos <- unmake_the_bed(file = "./hglft_female_cpg_38pos.bed")
male_cpg_38pos <- unmake_the_bed(file = "./hglft_male_cpg_38pos.bed")


wholesig_38 <- whole_cpg_38pos %>%
  left_join(supplementary2_Whole, by = c("CpG" = "probes"))
femalesig_38 <- female_cpg_38pos %>%
  left_join(supplementary2_Female, by = c("CpG" = "probes"))
malesig_38 <- male_cpg_38pos %>%
  left_join(supplementary2_Male, by = c("CpG" = "probes"))

#Create genomic ranges
library(GenomicRanges)
sig_grange <- function(df) {
  GRanges(
    seqnames = df$chr,
    ranges = IRanges(start = df$pos, end = df$pos),
    og_model_adj_P = df$adj.P.Val,
    deltaB = df$deltaB,
    CpG_ID = df$CpG
  )
}
wholesig38_grange <- sig_grange(wholesig_38)
malesig38_grange <- sig_grange(malesig_38)
femalesig38_grange <- sig_grange(femalesig_38)

locations_granges <- GRanges(
  seqnames = locations_df$seqnames,
  ranges = IRanges(start = locations_df$start, end = locations_df$end),
  gene = locations_df$symbol,
  type = locations_df$type,
  id = locations_df$id
)

#Reduce locations regions with overlap
gene_type <- split(locations_granges, paste(locations_granges$gene, locations_granges$type, sep = "_"))
reduced_locations_granges <- GenomicRanges::reduce(gene_type)
reduced_locations_gr <- unlist(reduced_locations_granges)
group_ids <- rep(names(reduced_locations_granges), elementNROWS(reduced_locations_granges))
mcols(reduced_locations_gr)$group_name <- group_ids

#overlap function
overlap_df <- function(data = df, locations = locations) {
  reduced_hits <- findOverlaps(data, locations)
  print(paste("Number of unique CpG hits:", length(unique(queryHits(reduced_hits)))))
  df2 <- data.frame(
    chr = as.character(seqnames(data)[queryHits(reduced_hits)]),
    position = start(data)[queryHits(reduced_hits)],
    og_model_adj_P = mcols(data)$og_model_adj_P[queryHits(reduced_hits)],
    deltaB = mcols(data)$deltaB[queryHits(reduced_hits)],
    region_chr = as.character(seqnames(locations)[subjectHits(reduced_hits)]),
    region_start = start(locations)[subjectHits(reduced_hits)],
    region_end = end(locations)[subjectHits(reduced_hits)],
    type = mcols(locations)$group_name[subjectHits(reduced_hits)],
    CpG_ID = mcols(data)$CpG_ID[queryHits(reduced_hits)]
    )
  df2 <- df2 %>%
    separate(
      col = type,
      into = c("gene", "region_type"),
      sep = "_hg38_genes_"
    )
}

reduced_overlap_whole <- overlap_df(data = wholesig38_grange, locations = reduced_locations_gr)
reduced_overlap_female <- overlap_df(data = femalesig38_grange, locations = reduced_locations_gr)
reduced_overlap_male <- overlap_df(data = malesig38_grange, locations = reduced_locations_gr)

#Function for Simplifying table for readability
write_readable_annotation <- function(overlap_df = overlap_data_frame) {
  region_labels <- c("1to5kb", "promoters", "5UTRs", "introns", "exons", "3UTRs")
  cpg_region_matrix <- overlap_df %>%
  dplyr::select(CpG_ID, gene, chr, position, og_model_adj_P, deltaB, region_type) %>%
  distinct() %>%  # Remove any duplicated overlaps
  mutate(overlap = TRUE) %>%
  pivot_wider(
    names_from = region_type,
    values_from = overlap,
    values_fill = FALSE  # Fill non-overlapping regions with FALSE
  ) %>% 
  mutate(adj_P_val = 
  og_model_adj_P
  ) %>%
  mutate(
  region_overlap = pmap_chr(
    across(any_of(region_labels)),  # select only existing columns
    ~ {
      labels <- region_labels[seq_along(list(...))]
      present <- unlist(list(...))
      paste(labels[which(present)], collapse = ",")
    }
  )
  ) %>% 
  mutate(
    region_overlap = gsub(",NA", "", region_overlap),
    region_overlap = gsub("NA,", "", region_overlap),
    region_overlap = gsub("NA", "", region_overlap),
    region_overlap = if_else(region_overlap == "", "none", region_overlap)
  ) %>%
  mutate(meth_status = case_when(
    deltaB < 0 & deltaB > -0.05 ~ "less DNAm",
    deltaB <= -0.05 ~ "Dif Hypomethylated",
    deltaB > 0 & deltaB < 0.05 ~ "more DNAm",
    deltaB >= 0.05 ~ "Dif Hypermethylated"
  )) %>%
  mutate (
    deltaB_rounded = round(deltaB, 3)
    ) %>%
  mutate(
    chr_num = as.numeric(gsub("chr", "", chr)) 
  ) %>%
  arrange(chr_num, position
  ) %>%
  transmute(
    CpG_ID,
    gene,
    chr,
    position,
    region_overlap,
    adj_P_val,
    meth_status,
    deltaB_rounded,
    )
  cpg_region_matrix <- as.data.frame(cpg_region_matrix)
}


whole_readable_annotation <- write_readable_annotation(reduced_overlap_whole)
write.csv(whole_readable_annotation, "./whole_readable_annotation.csv")
female_readable_annotation <- write_readable_annotation(reduced_overlap_female)
write.csv(female_readable_annotation, "./female_readable_annotation.csv")
male_readable_annotation <- write_readable_annotation(reduced_overlap_male)
write.csv(male_readable_annotation, "./male_readable_annotation.csv")


#Checking the overlap with the price anno 
price_overlap <- function(data = readable_annotation) {
  price_anno <- read.delim("/workspace/lab/wilsonslab/eyerk/2025_Lipid_Candidate_GSE_Info/GSE_annotations/GPL16304-47833_no_legend.tsv",header=TRUE)
  price_anno_simple <- price_anno %>%
    dplyr::select(ID, Closest_TSS_gene_name, Distance_closest_TSS)
  region_vs_price <- merge(data, price_anno_simple, by.x = "CpG_ID", by.y = "ID")
  region_vs_price <- region_vs_price %>%
    mutate(overlap = case_when(
      gene == Closest_TSS_gene_name ~ "Same", 
      gene != Closest_TSS_gene_name ~ " "
    )) %>%
    mutate(
      chr_num = as.numeric(gsub("chr", "", chr)) 
    ) %>%
    arrange(chr_num, position) 
}


whole_with_price <- price_overlap(whole_readable_annotation)
write.csv(whole_with_price, "./whole_region_vs_price38.csv")
female_with_price <- price_overlap(female_readable_annotation)
write.csv(female_with_price, "./female_region_vs_price38.csv")
male_with_price <- price_overlap(male_readable_annotation)
write.csv(male_with_price, "./male_region_vs_price38.csv")



















#Formatting into a nice table that shows this information in a nicer format 
cpg_region_matrix_df <- as.data.frame(cpg_region_matrix)
write.csv(cpg_region_matrix_df, "./cpg_region_matrix38.csv")







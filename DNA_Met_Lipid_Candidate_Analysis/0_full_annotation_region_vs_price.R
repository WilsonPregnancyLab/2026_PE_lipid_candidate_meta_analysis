#Adapted from mel-14's code Melanie_PhD_Code/IVF_MetaAnalysis/Significant_CpG_Annotation.R
# Map CpG sites from hg19 to hg38 and annotate to functional genomic regions (promoters, gene bodies, closest transcription start site)

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
library(rtracklayer)

# R 4.5.1 - done locally 
# Load hg38 annotations for promoters and gene bodies
setwd("/workspace/lab/wilsonslab/datalake-wilsonslab/2025_DNAm_Lipid_Candidate/Candidate_DNAm_Annotations")
annotations <- c("hg38_genes_promoters", "hg38_genes_5UTRs","hg38_genes_3UTRs","hg38_genes_introns","hg38_basicgenes")
locations <- build_annotations(genome = "hg38", annotations = annotations)
locations_df <- as.data.frame(locations)

gtf19_file <- "../../2025_RNA_Lipid_Candidate/genome_mapping/gencode.v37lift37.annotation.gtf"
gtf19_data <- import(gtf19_file)
gtf19_df <- as.data.frame(gtf19_data)

gtf_file <- "../../2025_RNA_Lipid_Candidate/genome_mapping/gencode.v48.primary_assembly.annotation.gtf"
gtf_data <- import(gtf_file)
gtf_df <- as.data.frame(gtf_data)

#R4.4.1 - done on Server 

price_anno <- read.delim("/workspace/lab/wilsonslab/eyerk/2025_Lipid_Candidate_GSE_Info/GSE_annotations/GPL16304-47833_no_legend.tsv",header=TRUE)

annotation19 <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
priceanno_19 <- as.data.frame(merge(price_anno, annotation19, by.x = "ID", by.y = "Name", all.x = TRUE))

#Function for making .bed files 
make_the_bed <- function(df = "df",
                         pos = "pos",
                         chr = "chr",
                         cpg = "cpg",
                         save_info = "/location/name.bed") {
  df2 <- df %>%
  mutate(posstart = as.numeric(pos) -1) %>%
  mutate(posend = paste(as.numeric(
    pos)
  )) %>%
  dplyr::select(chr, posstart, posend, cpg) %>%

  na.omit()

  write.table(df2,
    file = save_info,
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )
}


make_the_bed(priceanno_19, pos = "pos", chr = "chr", cpg = "ID", save_info = "./priceanno_19.bed")

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
#Convert using UCSC lift over (https://genome.ucsc.edu/cgi-bin/hgLiftOver)

#Load in the UCDC conversion to hg38 positions
price_cpg_38pos <- unmake_the_bed(file = "./hglft_price_cpg_38pos.bed")

price_38 <- price_cpg_38pos %>%
  left_join(priceanno_19, by = c("CpG" = "ID"))

#Create genomic ranges
library(GenomicRanges)
sig_grange <- function(df) {
  GRanges(
    seqnames = df$chr.x,
    ranges = IRanges(start = df$pos.x, end = df$pos.x),
    CpG_ID = df$CpG
  )
}
price38_grange <- sig_grange(price_38)

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
    #og_model_adj_P = mcols(data)$og_model_adj_P[queryHits(reduced_hits)],
    #deltaB = mcols(data)$deltaB[queryHits(reduced_hits)],
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

reduced_overlap_price <- overlap_df(data = price38_grange, locations = reduced_locations_gr)

#Function for Simplifying table for readability
write_readable_annotation <- function(overlap_df = overlap_data_frame) {
  region_labels <- c("1to5kb", "promoters", "5UTRs", "introns", "exons", "3UTRs")
  cpg_region_matrix <- overlap_df %>%
  dplyr::select(CpG_ID, gene, chr, position, region_type) %>%
  distinct() %>%  # Remove any duplicated overlaps
  mutate(overlap = TRUE) %>%
  pivot_wider(
    names_from = region_type,
    values_from = overlap,
    values_fill = FALSE  # Fill non-overlapping regions with FALSE
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
    )
  cpg_region_matrix <- as.data.frame(cpg_region_matrix)
}


price_readable_annotation <- write_readable_annotation(reduced_overlap_price)
write.csv(price_readable_annotation, "./price_readable_annotation.csv")



#need to merge Closest_TSS position from "price_cpg_38pos", match the gene_name from gtf_df, and the transcript_name from gtf_df
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
transcripts_38 <- transcripts(txdb)
tss_38 <- GenomicRanges::resize(transcripts_38, width = 1, fix = "start"
)
price_TSS_38_grange <- GRanges(
    seqnames = price_cpg_38pos$chr,
    ranges = IRanges(start = price_cpg_38pos$pos, end = price_cpg_38pos$pos),
    CpG_ID = price_cpg_38pos$CpG
  )

closest_TSS <- distanceToNearest(price_TSS_38_grange, tss_38)

closest_TSS_df <- data.frame(
  # Data from your original CpG query
  CpG_ID = mcols(price_TSS_38_grange)$CpG_ID[queryHits(closest_TSS)],
  CpG_Chr = as.character(seqnames(price_TSS_38_grange))[queryHits(closest_TSS)],
  CpG_Pos = start(price_TSS_38_grange)[queryHits(closest_TSS)],
  
  # Data extracted from the closest TSS subject
  Closest_TSS_Chr = as.character(seqnames(tss_38))[subjectHits(closest_TSS)],
  Closest_TSS_Pos = start(tss_38)[subjectHits(closest_TSS)],
  Closest_TSS_Transcript = mcols(tss_38)$tx_name[subjectHits(closest_TSS)],
  
  # The distance (in bp) calculated by the function
  Distance_Closest_TSS = mcols(closest_TSS)$distance
)

gtf_df_unique <- gtf_df[!duplicated(gtf_df$transcript_id), ]
closest_TSS_df <- merge(closest_TSS_df[,c("Closest_TSS_Transcript", "CpG_ID", "CpG_Chr", "CpG_Pos", "Closest_TSS_Chr", "Closest_TSS_Pos", "Distance_Closest_TSS")], gtf_df_unique[,c("transcript_id","gene_name", "gene_id")], by.x = "Closest_TSS_Transcript", by.y = "transcript_id")
colnames(closest_TSS_df)[colnames(closest_TSS_df) == "gene_name"] <- "Closest_TSS_gene_name"

#Check overlap with price anno
region_vs_price38 <- merge(price_readable_annotation, closest_TSS_df, by = "CpG_ID")
region_vs_price38 <- region_vs_price38 %>%
    mutate(overlap = case_when(
      gene == Closest_TSS_gene_name ~ "Same", 
      gene != Closest_TSS_gene_name ~ " "
    ))

write.csv(region_vs_price38, "./updated_full_annotation_region_vs_price38.csv")

# #Checking the overlap with the price anno 
# priceanno_overlap <- function(data = readable_annotation) {
#   price_anno <- read.delim("/workspace/lab/wilsonslab/eyerk/2025_Lipid_Candidate_GSE_Info/GSE_annotations/GPL16304-47833_no_legend.tsv",header=TRUE)
#   price_anno_simple <- price_anno %>%
#     dplyr::select(ID, Closest_TSS_gene_name, Distance_closest_TSS)
#   region_vs_price <- merge(data, price_anno_simple, by.x = "CpG_ID", by.y = "ID")
#   region_vs_price <- region_vs_price %>%
#     mutate(overlap = case_when(
#       gene == Closest_TSS_gene_name ~ "Same", 
#       gene != Closest_TSS_gene_name ~ " "
#     )) %>%
#     mutate(
#       chr_num = as.numeric(gsub("chr", "", chr)) 
#     ) %>%
#     arrange(chr_num, position) 
# }
# price_with_price <- priceanno_overlap(price_readable_annotation)

# write.csv(price_with_price, "./full_annotation_region_vs_price38.csv")





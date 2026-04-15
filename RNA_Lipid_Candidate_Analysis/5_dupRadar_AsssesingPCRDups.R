
#analyzeDuprates (in R)
library(dupRadar)

bam_files <- list.files(path = "/workspace/lab/wilsonslab/datalake-wilsonslab/2025_RNA_Lipid_Candidate/genome_mapping/markeddup_BAMs", pattern = ".bam", all.files = FALSE, full.names = TRUE)
gtf <- "/workspace/lab/wilsonslab/datalake-wilsonslab/2025_RNA_Lipid_Candidate/genome_mapping/gencode.v48.primary_assembly.annotation.gtf"

all_bam_dups <- list()

for (file in bam_files) {
    bamDups <-analyzeDuprates(file,gtf,stranded = 2, paired = TRUE, threads = 10, verbose = TRUE)
    srr <- sub("_.*", "", basename(file))
    output_filename <- paste0("analyzeDuprates_", srr, ".csv")
    write.csv(bamDups, output_filename)
    all_bam_dups[[output_filename]] <- bamDups
}

#retrieve y-intercept from regression to determine baseline duplication, using duprateExpFit
dupRadar_results <- data.frame()

for (sample_name in names(all_bam_dups)) {
    dupMatrix <- all_bam_dups[[sample_name]]
    dupFit <- duprateExpFit(dupMatrix)
    intercept_dup <- dupFit$intercept
    slope_dup <- dupFit$slope
    srr_csv <- sub("analyzeDuprates_", "", sample_name)
    srr_name <- sub(".csv", "", srr_csv)
    dupRadar_df <- data.frame(SRR_Number = srr_name, Intercept = intercept_dup, Slope = slope_dup)
    dupRadar_results <- rbind(dupRadar_results, dupRadar_df)
}

write.csv(dupRadar_results, "dupRadar_results.csv")
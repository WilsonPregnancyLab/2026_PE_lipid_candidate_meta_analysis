 ## dupRadar

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("dupRadar")
library(dupRadar)

#test AddOrReplaceReadGroups
java -jar picard.jar AddOrReplaceReadGroups I=trimmed_SRR10916911_Aligned.sortedByCoord.out.bam O=trimmed_SRR10916911_rg.bam RGID=SRR10916911 RGSM=SRR10916911 RGPL=Illumina RGLB=SRR10916911 RGPU=1
bamDuprm1a <- markDuplicates(dupremover="picard", bam="/workspace/lab/wilsonslab/datalake-wilsonslab/2025_RNA_Lipid_Candidate/genome_mapping/trimmed_SRR10916911_rg.bam", path="/workspace/lab/wilsonslab/eyerk/programs/picard/", rminput = FALSE)

BAM_FILES=("/workspace/lab/wilsonslab/datalake-wilsonslab/2025_RNA_Lipid_Candidate/genome_mapping/oldmapped_BAMs/"*.bam)

base=$(basename trimmed_SRR10916913_Aligned.sortedByCoord.out.bam .bam)
sample_name=${base:8:11}
java -jar picard.jar AddOrReplaceReadGroups I=/workspace/lab/wilsonslab/datalake-wilsonslab/2025_RNA_Lipid_Candidate/genome_mapping/oldmapped_BAMs/trimmed_SRR10916913_Aligned.sortedByCoord.out.bam O=./${sample_name}_rg.bam RGID=${sample_name} RGSM=${sample_name} RGPL=Illumina RGLB=${sample_name} RGPU=${sample_name}


#test analyzeDuprates
gtf <- "/workspace/lab/wilsonslab/datalake-wilsonslab/2025_RNA_Lipid_Candidate/genome_mapping/gencode.v48.primary_assembly.annotation.gtf"
file1 <- "/workspace/lab/wilsonslab/datalake-wilsonslab/2025_RNA_Lipid_Candidate/genome_mapping/markedduplicates_BAMs/trimmed_SRR10916911_Aligned.sortedByCoord.out.bam"
bamDup1 <-analyzeDuprates(file1,gtf,stranded = 2, paired = TRUE, threads = 4, verbose = TRUE, isDuplicate = NULL)

file2  <- "/workspace/lab/wilsonslab/datalake-wilsonslab/2025_RNA_Lipid_Candidate/genome_mapping/markedduplicates_BAMs/trimmed_SRR11498080_Aligned.sortedByCoord.out.bam"
bamDup2 <-analyzeDuprates(file2,gtf,stranded = 2, paired = TRUE, threads = 4, verbose = TRUE)

file3 <- "/workspace/lab/wilsonslab/datalake-wilsonslab/2025_RNA_Lipid_Candidate/genome_mapping/raw_BAMs/trimmed_SRR10916911_Aligned.sortedByCoord.out.bam"
bamDup3 <-analyzeDuprates(file3,gtf,stranded = 2, paired = TRUE, threads = 4, verbose = TRUE)

bamDup4 <- analyzeDuprates( bam = "/workspace/lab/wilsonslab/datalake-wilsonslab/2025_RNA_Lipid_Candidate/genome_mapping/raw_BAMs/trimmed_SRR10916911_Aligned.sortedByCoord.out.bam", gtf = "/workspace/lab/wilsonslab/datalake-wilsonslab/2025_RNA_Lipid_Candidate/genome_mapping/gencode.v48.primary_assembly.annotation.gtf", paired = TRUE,  stranded = 2, threads = 4)

file5 <- "/workspace/lab/wilsonslab/datalake-wilsonslab/2025_RNA_Lipid_Candidate/genome_mapping/trimmed_SRR10916911_rg_duprm.bam"
bamDup5 <-analyzeDuprates(file5,gtf,stranded = 2, paired = TRUE, threads = 4, verbose = TRUE)

#test plot
bamDup5Plot <- duprateExpDensPlot(DupMat=bamDup5) 
png("./bamDup5Plot.png", height = 9, width = 5, units = "in", res = 300)
par(mfrow=c(1,2))
bamDup5Plot <- duprateExpDensPlot(DupMat=bamDup5) 
dev.off()



#all bam files in directory
$BAM_FILES=("/workspace/lab/wilsonslab/datalake-wilsonslab/2025_RNA_Lipid_Candidate/genome_mapping/oldmapped_BAMs/"*.bam)
#AddOrReplaceReadGroups (required for MarkDuplicates) (in .sh)
for file in ${BAM_FILES[@]}; do
base=$(basename "$file" .bam)
sample_name=${base:8:11}
java -jar picard.jar AddOrReplaceReadGroups I=/workspace/lab/wilsonslab/datalake-wilsonslab/2025_RNA_Lipid_Candidate/genome_mapping/oldmapped_BAMs/trimmed_SRR10916913_Aligned.sortedByCoord.out.bam O=./${sample_name}_rg.bam RGID=${sample_name} RGSM=${sample_name} RGPL=Illumina RGLB=${sample_name} RGPU=${sample_name}


parallel -j 8 "basename=$(basename {} .bam) samplename=${base:8:11} java -jar picard.jar AddOrReplaceReadGroups I={} O=./RG_BAMs/${samplename}_rg.bam RGID=${samplename} RGSM=${samplename} RGPL=Illumina RGLB=${samplename} RGPU=${samplename}" ::: "${BAM_FILES[@]}"

#MarkDuplicates (in .sh)
RG_FILES=("/workspace/lab/wilsonslab/datalake-wilsonslab/2025_RNA_Lipid_Candidate/genome_mapping/markeddup_BAMs/RG_BAMs/"*.bam)
parallel -j 8 "base=$(basename {} .bam) sample_name=${base:8:11} java --enable-native-access=ALL-UNNAMED -jar picard.jar MarkDuplicates I={} O=./${sample_name}_dedup.bam M=./${sample_name}_dedup_metrics.txt TAGGING_POLICY=All" ::: "${RG_FILES[@]}"

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
#!/bin/bash


#Step 1: Install STAR

##Go to programs directory
mkdir /workspace/lab/wilsonslab/eyerk/STAR
cd /workspace/lab/wilsonslab/eyerk/STAR

##Install directory is the untarred folder after you download

##Download the STAR
wget https://github.com/alexdobin/STAR/archive/refs/tags/2.7.11b.tar.gz

##Extract Program
tar -zxvf 2.7.11b.tar.gz

##Move into extracted folder
cd STAR-2.7.11b

##Configure into installation directory (if binary is not already available)
cd STAR/source
make STAR

##Add the 'binaries (bin)' folder to yout PATH in .bashrc to allow command to run in any directory
echo 'export PATH="/workspace/lab/wilsonslab/eyerk/programs/STAR/STAR-2.7.11b/source:$PATH"' >> ~/.bashrc
##Reload bash
source ~/.bashrc

#Test it's working
which STAR
STAR --version



#Step 2: Generate Genome Indexes (help with mapping)

mkdir /workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/genome_mapping/
cd /workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/genome_mapping/

##Download reference genome sequences (FASTA files) and annotations (GTF file)

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/GRCh38.primary_assembly.genome.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.primary_assembly.annotation.gtf.gz
gunzip *.gz

##Generate Genome Index
mkdir ./genomeInd
STAR --runMode genomeGenerate --genomeDir ./genomeInd --genomeFastaFiles ./GRCh38.primary_assembly.genome.fa --sjdbGTFfile ./gencode.v48.primary_assembly.annotation.gtf --runThreadN 4 --sjdbOverhang 100 --sjdbGTFtagExonParentGene gene_id

#Step 3: Mapping Reads to Genome

##Create new folder to save mapped files
mkdir ./mapped_reads/mapped_PE/
TRIM_DIR="/workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/quality_control/FastP_Trimming"

PE_trimmed_1=("$TRIM_DIR"/*_1.fastq)
PE_trimmed_2=("$TRIM_DIR"/*_2.fastq)

parallel -j 4 --link ' base=$(basename {1} _1.fastq); STAR --runThreadN 4 --genomeDir ./genomeInd --readFilesIn {1} {2} --sjdbGTFfile ./gencode.v48.primary_assembly.annotation.gtf --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ./mapped_reads/mapped_PE/${base}_ ' ::: "${PE_trimmed_1[@]}" ::: "${PE_trimmed_2[@]}"


##Create new folder to transfer just the BAM files
mkdir /workspace/lab/wilsonslab/datalake-wilsonslab/2025_Lipid_PE/Diff_RNA_Expression/not_deduplicated_BAMS
mkdir /workspace/lab/wilsonslab/datalake-wilsonslab/2025_Lipid_PE/Diff_RNA_Expression/deduplicated_BAMS
mv ./mapped_reads/mapped_PE/*.bam /workspace/lab/wilsonslab/datalake-wilsonslab/2025_Lipid_PE/Diff_RNA_Expression/BAM_files/


#Step 4: Deduplicate BAM files to assess and remove PCR artefact

##Go to programs directory
mkdir /workspace/lab/wilsonslab/eyerk/picard
cd /workspace/lab/wilsonslab/eyerk/picard
wget https://github.com/broadinstitute/picard/releases/download/3.4.0/picard.jar
java -version
java -jar /workspace/lab/wilsonslab/eyerk/picard/picard.jar -h

cd /workspace/lab/wilsonslab/datalake-wilsonslab/2025_Lipid_PE/Diff_RNA_Expression/deduplicated_BAMS
BAM_FILES=("/workspace/lab/wilsonslab/datalake-wilsonslab/2025_Lipid_PE/Diff_RNA_Expression/BAM_files/"*.bam)
 
##EstimateLibraryComplexity
mkdir ./estlibcomplexity
parallel -j 4 "java -jar picard.jar EstimateLibraryComplexity I={} O=./estlibcomplexity/{.}_est_lib_complexity.csv" ::: "${BAM_FILES[@]}"

##MarkDuplicates (locates and tags duplciate reads in BAM file)
###output is new BAM file where duplicates are identified in SAM flags field for each read with hexadecimal value of 0x0400 (decimal value of 1024)
###can identify type of duplicate (optical vs PCR) using DT tag

parallel -j 4 "java -jar picard.jar MarkDuplicates I={1} O={.}_dedup.bam M={.}_dedup_metrics.csv --TAGGING_POLICY=All" ::: "${BAM_FILES[@]}"
















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
STAR --runMode genomeGenerate --genomeDir ./genomeInd --genomeFastaFiles ./GRCh38.primary_assembly.genome.fa --sjdbGTFfile ./gencode.v48.primary_assembly.annotation.gtf --runThreadN 4



#Step 3: Mapping Reads to Genome

##Create new folder to save mapped files
mkdir /path_for_mapped_files/

for file in /path_of_trimmed_fastqs/*.fastq; do
STAR --runThreadN 4 --genomeDir ./genomeInd --readFilesIn ${file} --outSAMtype BAM SortedByCoordinate --outFileNamePrefix /path_for_mapped_files/${file}; done

##Create new folder to transfer just the BAM files
mkdir /path_for_bam_files/
mv /path_for_mappped_files/*.bam /path_for_bam_files/
















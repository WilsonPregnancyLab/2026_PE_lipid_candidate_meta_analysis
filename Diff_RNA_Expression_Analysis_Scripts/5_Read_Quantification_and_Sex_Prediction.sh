#!/bin/bash


#Step 1: Install featureCounts

##Go to programs directory
mkdir /workspace/lab/wilsonslab/eyerk/featureCounts
cd /workspace/lab/wilsonslab/eyerk/featureCounts

##Install directory is the untarred folder after you download

##Download the featureCounts
wget https://sourceforge.net/projects/subread/files/subread-2.1.1/subread-2.1.1-source.tar.gz/

##Extract Program
tar -zxvf subread-2.1.1.tar.gz

##Move into extracted folder
cd subread-2.1.1

##Add the 'binaries (bin)' folder to yout PATH in .bashrc to allow command to run in any directory
echo 'export PATH="/workspace/lab/wilsonslab/eyerk/programs/subread/subread-2.1.1-Linux-x86_64/bin:$PATH"' >> ~/.bashrc
##Reload bash
source ~/.bashrc


#Step 2: Count Reads Using featureCounts

##need GTF annotation files and BAM files 

cd /workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/genome_mapping/

mkdir ./read_counts 

featureCounts -a ./gencode.v48.primary_assembly.annotation.gtf -p --countReadPairs -o ./read_counts.csv -T 4 ./mapped_BAMs/PE_mapped_BAMs/*.bam











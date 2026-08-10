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

#Move .txt files, and excluded .bam files (due to high duplciation) into new directory
cd /workspace/lab/wilsonslab/datalake-wilsonslab/2025_RNA_Lipid_Candidate/genome_mapping/markeddup_BAMs/
mkdir txt_files
mkdir excluded_bams

mv ./*.txt ./txt_files

# moved excluded .bam files into excluded_bams directory
mv ./SRR16505233_markdup.bam ./excluded_bams
mv ./SRR31033386_markdup.bam ./excluded_bams
mv ./SRR27882535_markdup.bam ./excluded_bams
mv ./SRR27882534_markdup.bam ./excluded_bams
mv ./SRR27882528_markdup.bam ./excluded_bams
mv ./SRR27882532_markdup.bam ./excluded_bams
mv ./SRR27882531_markdup.bam ./excluded_bams
mv ./SRR27882529_markdup.bam ./excluded_bams
mv ./SRR11498080_markdup.bam ./excluded_bams
mv ./SRR27882533_markdup.bam ./excluded_bams
mv ./SRR11498081_markdup.bam ./excluded_bams
mv ./SRR27882527_markdup.bam ./excluded_bams
mv ./SRR27882530_markdup.bam ./excluded_bams
mv ./SRR27882536_markdup.bam ./excluded_bams
mv ./SRR24900370_markdup.bam ./excluded_bams
mv ./SRR31033382_markdup.bam ./excluded_bams
mv ./SRR31033384_markdup.bam ./excluded_bams
mv ./SRR31033383_markdup.bam ./excluded_bams
mv ./SRR31033385_markdup.bam ./excluded_bams
mv ./SRR31033376_markdup.bam ./excluded_bams
mv ./SRR31033387_markdup.bam ./excluded_bams


cd /workspace/lab/wilsonslab/datalake-wilsonslab/2025_RNA_Lipid_Candidate/

mkdir ./read_counts 

featureCounts -a ../genome_mapping/gencode.v48.primary_assembly.annotation.gtf -p --countReadPairs -s 2 -o ./read_counts.txt -T 10 --verbose ../genome_mapping/markeddup_BAMs/*.bam

#Sex Prediction using SexChrLab
git clone https://github.com/SexChrLab/SexInference.git
cd ../
mkdir sex_prediction
cd sex_prediction

#need to double-check that txtimport.R have transcript ID and gene ID that match "gencode.v48.primary_assembly.annotation.gtf"
#the software uses gencode.v29.annotation.gtf by default (lines 38-43 of txtimport.R)
#Rscript tximport_v48.R ../genome_mapping/gencode.v48.primary_assembly.annotation.gtf "./sample_name_SexInference.csv" "../read_counts.txt" "./data_for_regression_2025.csv"
#Rscript sex_inference_model.R "./data_for_regression_2025.csv"
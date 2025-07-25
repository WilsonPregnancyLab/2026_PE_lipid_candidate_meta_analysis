!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

###FAST_QC


## First, we will run FastQC which helps identify whether adaptor sequences are present in RNA-seq data, low-quality reads/base pairing, so that we can trim them
## We will run FastQC, then trim, then run FastQC again to ensure the trimming worked

#identify input and output directories we want to save quality control files (.HTML) (here, based on platform type)
fastq_sets_1=(
  "/workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/raw_data/Illumina_NovaSeq_6000/GSE255126 /workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/quality_control/Illumina_NovaSeq_6000/GSE255126/"
  "/workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/raw_data/Illumina_NovaSeq_6000/GSE279757 /workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/quality_control/Illumina_NovaSeq_6000/GSE279757/")

fastq_sets_2=(
  "/workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/raw_data/Illumina_NovaSeq_6000/GSE186257 /workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/quality_control/Illumina_NovaSeq_6000/GSE186257/"
  "/workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/raw_data/Illumina_NovaSeq_6000/GSE234729 /workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/quality_control/Illumina_NovaSeq_6000/GSE234729/"
  "/workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/raw_data/Illumina_NovaSeq_6000/GSE218039 /workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/quality_control/Illumina_NovaSeq_6000/GSE218039/")

fastq_sets_3=(
  "/workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/raw_data/Illumina_NextSeq_2000/GSE204835 /workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/quality_control/Illumina_NextSeq_2000/GSE204835/"

  "/workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/raw_data/Illumina_HiSeq_2000/GSE114691 /workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/quality_control/Illumina_HiSeq_2000/GSE114691/"

  "/workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/raw_data/Illumina_HiSeq_4000/GSE143953 /workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/quality_control/Illumina_HiSeq_4000/GSE143953/")

fastq_sets_4=(
  "/workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/raw_data/Illumina_HiSeq_2500/GSE148241 /workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/quality_control/Illumina_HiSeq_2500/GSE148241/"
  "/workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/raw_data/Illumina_HiSeq_2500/GSE203507 /workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/quality_control/Illumina_HiSeq_2500/GSE203507/"
)

#Automate FASTQ file quality control
## the function was run in 4 separate windows for each of the fastq_sets to reduce the time to run

for pair in "${fastq_sets_1[@]}"; do
  # Split input and output
  input_dir=$(echo "$pair" | cut -d' ' -f1)
  output_dir=$(echo "$pair" | cut -d' ' -f2)

  for fastq in "$input_dir"/*.fastq*; do
    echo "  > Running FASTQC on: $fastq"
    fastqc --outdir "$output_dir" -t 4 "$fastq"
    echo "> Finished $fastq"
  done
  echo "Finished folder: $input_dir"
done
echo " All datasets complete."




###FASTQ_SCREEN

## Next, will run FastQScreen which ensures the RNAseq data maps to the human genome and not other species (bacteria, mouse, etc.)
## Could integrate this with previous code

#Download FastQScreen here: https://github.com/StevenWingett/FastQ-Screen/archive/refs/tags/v0.16.0.tar.gz

tar -xzf FastQ-Screen-0.16.0.tar.gz
wget https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/genome_locations.txt
$RUN ./fastq_screen --get_genomes

# Test 
# $RUN ./fastq_screen --conf FastQ_Screen_Genomes/fastq_screen.conf \
           /workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/raw_data/Illumina_NovaSeq_6000/GSE255126/SRR27882527_1.fastq \
           --outdir /workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/quality_control/FastQ_Screen/Illumina_NovaSeq_6000/GSE255126/

#identify input and output directories we want to save FASTQ Screen (.HTML) (here, based on platform type)
fastq_sets_1=(
  "/workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/raw_data/Illumina_NovaSeq_6000/GSE255126 /workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/quality_control/FastQ_Screen/Illumina_NovaSeq_6000/GSE255126/"
  "/workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/raw_data/Illumina_NovaSeq_6000/GSE279757 /workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/quality_control/FastQ_Screen/Illumina_NovaSeq_6000/GSE279757/")

fastq_sets_2=(
  "/workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/raw_data/Illumina_NovaSeq_6000/GSE186257 /workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/quality_control/FastQ_Screen/Illumina_NovaSeq_6000/GSE186257/"
  "/workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/raw_data/Illumina_NovaSeq_6000/GSE234729 /workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/quality_control/FastQ_Screen/Illumina_NovaSeq_6000/GSE234729/"
  "/workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/raw_data/Illumina_NovaSeq_6000/GSE218039 /workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/quality_control/FastQ_Screen/Illumina_NovaSeq_6000/GSE218039/")

fastq_sets_3=(
  "/workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/raw_data/Illumina_NextSeq_2000/GSE204835 /workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/quality_control/FastQ_Screen/Illumina_NextSeq_2000/GSE204835/"

  "/workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/raw_data/Illumina_HiSeq_2000/GSE114691 /workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/quality_control/FastQ_Screen/Illumina_HiSeq_2000/GSE114691/"

  "/workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/raw_data/Illumina_HiSeq_4000/GSE143953 /workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/quality_control/FastQ_Screen/Illumina_HiSeq_4000/GSE143953/")

fastq_sets_4=(
  "/workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/raw_data/Illumina_HiSeq_2500/GSE148241 /workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/quality_control/FastQ_Screen/Illumina_HiSeq_2500/GSE148241/"
  "/workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/raw_data/Illumina_HiSeq_2500/GSE203507 /workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/quality_control/FastQ_Screen/Illumina_HiSeq_2500/GSE203507/"
)


#Run FASTQ_Screen
CONFIG_PATH="/workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/quality_control/FastQ_Screen/FastQ-Screen-0.16.0/FastQ_Screen_Genomes/fastq_screen.conf"

for pair in "${fastq_sets_4[@]}"; do
  # Split input and output
  input_dir=$(echo "$pair" | cut -d' ' -f1)
  output_dir=$(echo "$pair" | cut -d' ' -f2)

  for fastq in "$input_dir"/*.fastq*; do
    echo "  > Running FASTQ Screen on: $fastq"
    ./fastq_screen --conf "$CONFIG_PATH" --outdir "$output_dir" "$fastq"
    echo "> Finished FASTQ Screen on $fastq"
  done
  echo "Finished FASTQ Screen on folder: $input_dir"
done
echo " All datasets complete."





#FASTP - Adapter Trimming

##identify input and output directories we want to save FASTQ Screen (.HTML) (here, based on platform type)
PE_input=(
  "/workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/raw_data/Illumina_NovaSeq_6000/GSE186257"
  "/workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/raw_data/Illumina_NovaSeq_6000/GSE234729"
  "/workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/raw_data/Illumina_NovaSeq_6000/GSE218039"
  "/workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/raw_data/Illumina_HiSeq_2500/GSE148241"
  "/workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/raw_data/Illumina_HiSeq_4000/GSE143953")

###CHECK GSE186257 AFTER
"/workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/raw_data/Illumina_NovaSeq_6000/GSE255126"
  "/workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/raw_data/Illumina_NovaSeq_6000/GSE279757"


PE_fastqs_1=()
PE_fastqs_2=()

for dir in "${PE_input[@]}"; do
  PE_fastqs_1+=("$dir"/*_1.fastq)
  PE_fastqs_2+=("$dir"/*_2.fastq)
done

SE_input=(
  "/workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/raw_data/Illumina_HiSeq_2000/GSE114691"
  "/workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/raw_data/Illumina_HiSeq_2500/GSE203507")

SE_fastqs=()
for dir in "${SE_input[@]}"; do
  SE_fastqs+=("$dir"/*.fastq)
done

##Install fastp
wget http://opengene.org/fastp/fastp
chmod a+x ./fastp
echo 'export PATH="/workspace/lab/wilsonslab/eyerk/programs/fastp-0.23.4/bin:$PATH"' >> ~/.bashr

##Install parallel
wget http://ftpmirror.gnu.org/parallel/parallel-latest.tar.bz2
tar -xvjf parallel-latest.tar.bz2
cd parallel-20250722
./configure --prefix=/workspace/lab/wilsonslab/eyerk/programs/parallel-20250722/
make
make install
echo 'export PATH="/workspace/lab/wilsonslab/eyerk/programs/parallel-20250722/bin:$PATH"' >> ~/.bashrc

#simple usage, -i indicates first input file in paired-end, -I indicates second input file in paired-end
#fastp -i in.R1.fq.gz -I in.R2.fq.gz -o out.R1.fq.gz -O out.R2.fq.gz -h --failed_out /workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/quality_control/FastP_Trimming/Failed_Reads

##fastp- trimming the adapters and appending UMI to sample same
## -i input of read 1
## -I input of read 2
## -o output of read 1
## -O output of read 2
## --linksm read 1 and read 2
## fastp automatically gzips outputs

cd /workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/quality_control/FastP_Trimming/

###Samples with Single-Ends and need UMI removal (GSE204835) 
parallel -j 4 'fastp -i {} -o ./trimmed_{/} --umi --umi_loc=per_read --umi_len=6  -j ./{/.}.json -h ./{/.}.html' ::: /workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/raw_data/Illumina_NextSeq_2000/GSE204835/*.fastq

###Samples with Single-Ends 
parallel -j 4 "fastp -i {} -o ./trimmed_{/} -j ./{/.}.json -h ./{/.}.html" ::: "${SE_fastqs[@]}"

###Samples with Paired-Ends 
parallel -j 2 --link "fastp -i {1} -I {2} -o ./trimmed_{1/}  -O ./trimmed_{2/} -j ./{1/.}.{2/.}.json -h ./{1/.}.{2/.}.html" ::: ${PE_fastqs_1[@]} ::: ${PE_fastqs_2[@]} --tmpdir /workspace/lab/wilsonslab/eyerk/temp


#FASTQC - Post-trimming Quality Control

cd /workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/quality_control/FastP_Trimming/

fastqc --outdir /workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/quality_control/FastQC_aftertrim/ -t 4 *.fastq
echo " All datasets complete."



#MultiQC - Consolidating FastQC Results Before and After Trimming, FastScreen and Fastp

multiqc /workspace/lab/wilsonslab/eyerk/2025_RNA_Lipid_Candidate/quality_control/


















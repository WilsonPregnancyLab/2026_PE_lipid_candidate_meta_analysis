#Download SRA Toolkit onto server
cd /workspace/lab/wilsonslab/eyerk/
mkdir ./SRA_Toolkit/
mkdir ./temp/
cd ./SRA_Toolkit/
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.3.0/sratoolkit.3.3.0-alma_linux64.tar.gz #version 3.2.1
tar -xvzf sratoolkit.3.2.1-alma_linux64.tar.gz

#Configure SRA Toolkit
cd ./sratoolkit.3.2.1-alma_linux64/bin/
./vdb-config -i
## You will be taken into an interactive terminal. Use the 'tab' button to navigate to the "CACHE" tab.
## Change "location of user-repository" to the temp folder you created above. Save changes and exit terminal.

#Download Sequencing Data
## Go to SRA Run Selector website for the studies to be downloaded, download the "Accession List" from the "Select" box

#location of SRA toolkit workspace
SRA_TOOLKIT_PATH="/workspace/lab/wilsonslab/eyerk/SRA_Toolkit/sratoolkit.3.2.1-alma_linux64/bin"
$TEMP_DIR=""

#make directories to save FASTQ files (here, based on platform type)
FASTQ_DIR_NS6000_GSE255126="/workspace/lab/wilsonslab/eyerk/raw_data/PE_RNA_seq_data/Illumina_NovaSeq_6000/GSE186257/"
FASTQ_DIR_NS6000_GSE279757="/workspace/lab/wilsonslab/eyerk/raw_data/PE_RNA_seq_data/Illumina_NovaSeq_6000/GSE234729/"
FASTQ_DIR_NS6000_GSE186257="/workspace/lab/wilsonslab/eyerk/raw_data/PE_RNA_seq_data/Illumina_NovaSeq_6000/GSE186257/"
FASTQ_DIR_NS6000_GSE234729="/workspace/lab/wilsonslab/eyerk/raw_data/PE_RNA_seq_data/Illumina_NovaSeq_6000/GSE234729/"
FASTQ_DIR_HS2500_GSE148241="/workspace/lab/wilsonslab/eyerk/raw_data/PE_RNA_seq_data/Illumina_HiSeq_2500/GSE148241/"
FASTQ_DIR_HS4000_GSE143953="/workspace/lab/wilsonslab/eyerk/raw_data/PE_RNA_seq_data/Illumina_HiSeq_4000/GSE143953/"

#define the directory you saved the accession list and the corresponding directory the fastq files will be saved ("location/of/accession/list *space* $FASTQ_DIR")

accession_sets=(
  "/workspace/lab/wilsonslab/eyerk/raw_data/PE_RNA_seq_data/Illumina_NovaSeq_6000/SRR_Acc_List_GSE255126.txt $FASTQ_DIR_NS6000_GSE255126"
  "/workspace/lab/wilsonslab/eyerk/raw_data/PE_RNA_seq_data/Illumina_NovaSeq_6000/SRR_Acc_List_GSE279757.txt $FASTQ_DIR_NS6000_GSE279757"
  "/workspace/lab/wilsonslab/eyerk/raw_data/PE_RNA_seq_data/Illumina_NovaSeq_6000/SRR_Acc_List_GSE186257.txt $FASTQ_DIR_NS6000_GSE186257"
  "/workspace/lab/wilsonslab/eyerk/raw_data/PE_RNA_seq_data/Illumina_NovaSeq_6000/SRR_Acc_List_GSE234729.txt $FASTQ_DIR_NS6000_GSE234729"
  "/workspace/lab/wilsonslab/eyerk/raw_data/PE_RNA_seq_data/Illumina_HiSeq_4000/SRR_Acc_List_GSE143953.txt $FASTQ_DIR_HS4000_GSE143953"
  "/workspace/lab/wilsonslab/eyerk/raw_data/PE_RNA_seq_data/Illumina_HiSeq_2500/SRR_Acc_List_GSE148241.txt $FASTQ_DIR_HS2500_GSE148241"
)


#Automate accession file downloads

for entry in "${accession_sets[@]}"; do
  #Divide into accession file and output directory
  accession_file=$(echo "$entry" | cut -d' ' -f1)
  output_dir=$(echo "$entry" | cut -d' ' -f2)
  
  echo ">>> Processing $accession_file"

  while read -r accession; do 
    echo " >>> Processing $accession"
    "$SRA_TOOLKIT_PATH/prefetch" --output-directory "$output_dir" "$accession"
    # -t defines a directory where temporary files can be stored
    "$SRA_TOOLKIT_PATH/fasterq-dump" --split-files --outdir "$output_dir" -t /workspace/lab/wilsonslab/eyerk/raw_data/PE_RNA_seq_data/temp "$accession"
    echo ">>> Finished $accession"
  done < "$accession_file"
  echo ">>> Completed $accession_file"

done


#Mapping RNA-Seq Data to Transcriptome Only (as opposed to genome)

salmon --version #salmon 1.10.3

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.transcripts.fa.gz

salmon index -t gencode.v48.transcripts.fa.gz -i gencode_v48_transcripts_index --gencode

mkdir ./quants

for r1 in /workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/quality_control/multiqc_after_processing/FastP_Trimming/*_1.fastq
do
    base=$(basename "$r1" _1.fastq)
    r2="${r1/_1/_2}"
    echo "Processing sample ${base}"
    salmon quant -i gencode_v48_transcripts_index -l A \
         -1 "$r1" \
         -2 "$r2" \
         -p 12 --validateMappings -o quants/${base}_quant --numBootstraps 30
done 

gunzip /workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2026_Alternative_Splicing/gencode.v48.transcripts.fa.gz

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.annotation.gtf.gz
gunzip gencode.v48.annotation.gtf.gz

#Differential Alternative Splicing Analysis

wget http://rnaseq-mats.sourceforge.net/rMATS.4.0.2.tgz
tar -xzvf rmats_turbo_v4_3_0.tar.gz

##Install gsl
wget https://ftp.gnu.org/gnu/gsl/gsl-2.5.tar.gz
tar -xzf gsl-2.5.tar.gz
cd gsl-2.5
./configure --prefix=$HOME/programs/local
make -j$(nproc)
make install

echo 'export PATH="/workspace/lab/wilsonslab/eyerk/programs/local/bin:$PATH"' >> ~/.bashrc
echo 'export LD_LIBRARY_PATH="/workspace/lab/wilsonslab/eyerk/programs/local/lib:$LD_LIBRARY_PATH"' >> ~/.bashrc
source ~/.bashrc

#!/bin/bash
LOCAL_BIN="/workspace/lab/wilsonslab/eyerk/programs/local/bin"
LOCAL_LIB="/workspace/lab/wilsonslab/eyerk/programs/local/lib"
RMATS_DIR="/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2026_Alternative_Splicing/rmats_turbo_v4_3_0"

#set default python library as personal library as opposed to system library
export LD_LIBRARY_PATH="$LOCAL_LIB:$LD_LIBRARY_PATH"
export PYTHONPATH="${HOME}/.local/lib/python3.9/site-packages:$PYTHONPATH"

# 2. SYSTEM LIBRARY LINKING
mkdir -p "$LOCAL_LIB"
ln -sf /usr/lib64/liblapack.so.3 "$LOCAL_LIB/liblapack.so"
ln -sf /usr/lib64/libblas.so.3 "$LOCAL_LIB/libblas.so"

# install dependencies
python3 -m pip install --user Cython setuptools numpy

# fix python 3.12 compatibility
cd "$RMATS_DIR/rMATS_pipeline"
sed -i 's/from distutils.core import setup/from setuptools import setup, Extension/g' setup.py
sed -i '/from distutils.extension import Extension/d' setup.py
python3 setup.py build_ext --inplace

cd /workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2026_Alternative_Splicing/rmats_turbo_v4_3_0/
./build_rmats
cp rMATS_pipeline/rmats/rmatspipeline*.so .

python3 rmats.py --version

#Run from any directory
echo "alias rmats='python3 /workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2026_Alternative_Splicing/rmats_turbo_v4_3_0/rmats.py'" >> ~/.bashrc
source ~/.bashrc


#in R
metadata <- read.csv("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/metadata/RNA_Lipid_Candidate_Metadata.csv")
metadata <- metadata[metadata$exclude != "exclude",]
metadata$BAM_filename <- paste0("/workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/genome_mapping/mapped_reads/mapped_PE/trimmed_", metadata$Run, "_Aligned.sortedByCoord.out.bam")
PE_Bam <- metadata[metadata$disease_group == "PE", ]
PE_Bam_F <- metadata[metadata$disease_group == "PE" & metadata$predicted_fetal_sex == "F", ]
PE_Bam_M <- metadata[metadata$disease_group == "PE" & metadata$predicted_fetal_sex == "M", ]

cat(PE_Bam$BAM_filename, sep = ",", file = "PE_BAM_filename.txt")
cat(PE_Bam_F$BAM_filename, sep = ",", file = "PE_BAM_F_filename.txt")
cat(PE_Bam_M$BAM_filename, sep = ",", file = "PE_BAM_M_filename.txt")

Cont_Bam <- metadata[metadata$disease_group == "Control", ]
Cont_Bam_F <- metadata[metadata$disease_group == "Control" & metadata$predicted_fetal_sex == "F", ]
Cont_Bam_M <- metadata[metadata$disease_group == "Control" & metadata$predicted_fetal_sex == "M", ]

cat(Cont_Bam$BAM_filename, sep = ",", file = "CONt_BAM_filename.txt")
cat(Cont_Bam_F$BAM_filename, sep = ",", file = "CONt_BAM_F_filename.txt")
cat(Cont_Bam_M$BAM_filename, sep = ",", file = "CONt_BAM_M_filename.txt")

#rMATS Combined Fetal Sex
rmats --b1 /path/to/b1.txt --b2 /path/to/b2.txt --gtf /path/to/the.gtf -t paired --readLength 50 --nthread 12 --od /path/to/output --tmp /path/to/tmp_output

python3 rmats.py --b1 /workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2026_Alternative_Splicing/alt_splic_analysis/PE_BAM_filename.txt --b2 /workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2026_Alternative_Splicing/alt_splic_analysis/CONt_BAM_filename.txt --gtf /workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/genome_mapping/gencode.v48.primary_assembly.annotation.gtf -t paired --readLength 141 --variable-read-length --allow-clipping --nthread 12 --od /workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2026_Alternative_Splicing/alt_splic_analysis/output --tmp /workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2026_Alternative_Splicing/alt_splic_analysis/tmp
python3 rmats.py --b1 /workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2026_Alternative_Splicing/alt_splic_analysis/PE_BAM_F_filename.txt --b2 /workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2026_Alternative_Splicing/alt_splic_analysis/CONt_BAM_F_filename.txt --gtf /workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/genome_mapping/gencode.v48.primary_assembly.annotation.gtf -t paired --readLength 141 --variable-read-length --allow-clipping --nthread 12 --od /workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2026_Alternative_Splicing/alt_splic_analysis/output_F --tmp /workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2026_Alternative_Splicing/alt_splic_analysis/tmp_F
python3 rmats.py --b1 /workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2026_Alternative_Splicing/alt_splic_analysis/PE_BAM_M_filename.txt --b2 /workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2026_Alternative_Splicing/alt_splic_analysis/CONt_BAM_M_filename.txt --gtf /workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2025_RNA_Lipid_Candidate/genome_mapping/gencode.v48.primary_assembly.annotation.gtf -t paired --readLength 141 --variable-read-length --allow-clipping --nthread 12 --od /workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2026_Alternative_Splicing/alt_splic_analysis/output_M --tmp /workspace/lab/wilsonslab/eyerk/PE_Lipid_Meta-analysis/2026_Alternative_Splicing/alt_splic_analysis/tmp_M


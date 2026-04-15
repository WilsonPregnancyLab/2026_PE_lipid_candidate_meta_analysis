
#Step 4: Assess alignment quality

cd /workspace/lab/wilsonslab/eyerk/programs/
mkdir ./qualimap
cd ./qualimap
wget https://bitbucket.org/kokonech/qualimap/downloads/qualimap_v2.3.zip
unzip qualimap_v2.3.zip
echo 'export PATH="/workspace/lab/wilsonslab/eyerk/programs/qualimap/qualimap_v2.3:$PATH"' >> ~/.bashrc
source ~/.bashrc

cd /workspace/lab/wilsonslab/datalake-wilsonslab/2025_RNA_Lipid_Candidate/genome_mapping/alignment_qualimap/
qualimap multi-bamqc -d /workspace/lab/wilsonslab/datalake-wilsonslab/2025_RNA_Lipid_Candidate/genome_mapping/alignment_qualimap/qualimap_input_file.txt -outdir ./ -r

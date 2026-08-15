This readme file was generated on 2026-08-10 by Kriesha Eyer.

# General Information

Title of Dataset: 2026_PE_lipid_candidate_meta_analysis

Author/Principal Investigator Information
Name: Samantha L. Wilson
ORCID:
Institution: McMaster University 
Address: 1200 Main Street W, Hamilton, Ontario L8N3Z5
Email: wilsos66@mcmaster.ca

Author/Associate or Co-investigator Information
Name: Kriesha Eyer
ORCID:
Institution: McMaster University 
Address: 1200 Main Street W, Hamilton, Ontario L8N3Z5
Email: eyerk@mcmaster.ca

Author/Alternate Contact Information
Name: Melanie Lemaire
ORCID:
Institution: McMaster University 
Address: 1200 Main Street W, Hamilton, Ontario L8N3Z5
Email: lemairem@mcmaster.ca

Date of data collection: 2024-10-09 - 2026-06-01 

Geographic location of data collection: Hamilton, Ontario, Canada

# Sharing/Access Information

### Datasets Used for Lipid Candidate Gene Search: 
- GAD Gene-Disease Associations Dataset - https://maayanlab.cloud/Harmonizome/dataset/GAD+Gene-Disease+Associations
- GO Biological Process Annotations 2023 Dataset - https://maayanlab.cloud/Harmonizome/dataset/GO+Biological+Process+Annotations+2023
-  GO Cellular Component Annotations 2023 Dataset - https://maayanlab.cloud/Harmonizome/dataset/GO+Cellular+Component+Annotations+2023
- Reactome Pathways 2024 Dataset - https://maayanlab.cloud/Harmonizome/dataset/Reactome+Pathways+2024
  
### Datasets Downloaded from Gene Expression Omnibus: 
RNA-seq Data
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE255126
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE148241
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE143953
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186257
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE234729
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE279757

Affymetrix Human Gene 1.0 ST Array Data
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75010

DNA Methylation Data
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125605
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100197
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE98224
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75196

# Data & File Overview

All project scripts to run analyses and produce. Order of read me is the order scipts should be run. R version 4.6.0 was used to run all scripts.

* **Scripts/** (all execution scripts)
    * **1_RNA_Lipid_Candidate_Analysis/** 
        * `0_lipid_candidate_gene_extraction.R`: Identify lipid-related genes from Reactome, Gene-Disease Associations and lipid databases, and papers based on key lipid-related terms.
        * `1_SRA_NCBI_RNA_Seq_data_download.sh`: Download RNA-sequencing data from NCBI GEO using SRA Toolkit.
        * `2_metadata_extraction_RNAseq.R`: Extract clinical and technical data for all RNA-seq samples using GEOquery.

        * `3_QC_and_Trimming.sh`: Perform quality control and read trimming using FastQC and fastp.
        * `4_Mapping_Reads_to_Genome.sh`: Mapping RNA-seq reads to gencodev48 using STAR, and marking duplicate reads in BAM files to remove PCR artifact.
        * `5_dupRadar_AsssesingPCRDups.R`: Assessing duplicates using dupRadar.
        * `6_Fetal_Sex_Prediction.R`: Predicting fetal sex of samples using SexInference package from SexChrLab. Used this information as a covariate in model.
        * `7_Read_Quantification_featureCounts.sh`: Count number of reads mapping to each gene using featureCounts.
        * `8a_DESeq2_Lip_Candidate_differential_exp_analysis.R`: Use DESeq2 to perform differential lipid-candiate gene expression analysis between PE and control samples. Plot volcano plots using ggplot2.
        * `8b_DESeq2_WG_differential_exp_analysis.R`: Use DESeq2 to perform genome-wide differential gene expression analysis between PE and control samples. Plot volcano plots using ggplot2.
        * `9_PCA_plot_RNA_seq.R`: Construct PCA plots before and after normalization with DESeq2 and ggplot2.
        * `10a_Pathway_Enrichment_RNA_Seq.R`: Use ErmineJ to identify the most common pathways that have differential gene expression in lipid-candidate genes.
        * `10b_WG_Pathway_Enrichment_RNAseq.R`: Use ErmineJ to identify the most common pathways that have differential gene expression across genome.
        * `11_RNA_Cell_Deconvolution.R`: Use MuSiC to identify differences in cell-type proportion between PE and control samples based on RNA-seq data.
        * `12_Parsing_Out_Sex_Differences.R`: Identify genes in sex-stratified analyses unique to each sex, and reperform DESeq2 analysis while adding interaction term (fetal_sex*pathology_group) to identify sex-specific pathology.
        * `13_Salmon_and_Alternative_Splicing.sh`: Map RNA-seq reads to gencodev48, this time at transcript-level using Salmon. Install rMATS.
        * `14_Differential_Isoform_Analysis.R`: Perform differential transcript/isoform expression analysis using Swish. Run rMATS to identify sites of alternative splicing in PE compared to controls.
        * `15a_Isoform_Switch_Analysis_combinedfetalsex.R`: Perform isoform switch analysis in the combined fetal-sex population using IsoformSwitchAnalyzeR to identify whether reciprocal changes in isoform expression occured in PE. 
        * `15b_Isoform_Switch_Analysis_stratfetalsex.R`: Perform isoform switch analysis in sex-stratified populations using IsoformSwitchAnalyzeR to identify whether reciprocal changes in isoform expression occured in PE.
        * **Affymetrix_Validation/**
            * `1_Affymetrix_QC_Normalization.R`: Download affymetrix dataset and perform quality control and normalization using arrayQualityMetrics.
            * `2a_Affymetrix_Lip_PE_RNA_Diff_Exp.R`: Use delta exxpression values between control and preeclamptic samples to identify significant differentially expressed lipid-candidate genes. Plot using ggplot2.
            * `2b_Affymetrix_WG_RNA_Diff_Exp.R`: Use delta exxpression values between control and preeclamptic samples to identify significant differentially expressed genes genome-wide. Plot using ggplot2.

    * **2_DNA_Met_Lipid_Candidate_Analysis/**
        * `1_full_annotation_region_vs_price.R`: Use price annotation file and gtf file to map CpG probes to gene location.
        * `2_GEO_dataset_download.sh`: Download all DNAm datasets from NCBI GEO. 
        * `3_metadata_extraction_DNAm.R`: Extract clinical and technical data for all DNAm sample sets using GEOquery.
        * `4_fetal_sex_prediction.R`: Use minfi package to estimate fetal sex for all DNAm samples as a quality control measure and to include fetal sex as a covariate. 
        * `5_normalization_probefiltering.R`: Exclude samples from other diseases. Normalize data with adjustedFunnorm. Exclude probes with bad detection P-values, missing beta-values, SNP probes, cross-hybridizing probes and non-variable placental probes. 
        * `6_PCA_plot_DNAm.R`: Construct PCA plots before and after normalization with adjustedFunnorm using missMDA and ggplot2.
        * `7a_Diff_DNA_Met_Analysis.R`: Use delta beta values between control and preeclamptic samples to identify significant differentially methylated CpGs in lipid-related genes. Use M-values to perform linear modeling for visualization of results using volcano plots.
        * `7b_WG_Diff_DNA_Met_Analysis.R`: Use delta beta values between control and preeclamptic samples to identify significant differentially methylated CpGs across genome. Use M-values to perform linear modeling for visualization of results using volcano plots.
        * `8a_DMR_lipid_candidate_genes.R`: Use DMRcate to conduct differential regional DNAm analysis in lipid-candidate genes.
        * `8b_DMR_whole_genome_analysis.R`: Use DMRcate to conduct differential regional DNAm analysis across genome.
        * `9_pathway_enrichment_analysis.R`: Use ErmineJ to identify the most common pathways that have differential DNAm.
        * `10_DNAm_Cell_Deconvolution.R`: Use planet and EpiDish to identify differences in cell-type proportion between PE and control samples based on DNAm data.

    * **3_Combined_DNA_RNA_Analysis/**
        * `DNAm_RNA_affy_overlap.R`: Identify overlap between differential expression and DNAm across all platforms.
        * `Lip_DNAm_RNA_Upset_Plot.R`: Use UpSetR to create upset plot showing number of lipid-related genes with differential DNAm and expression across platforms.  
        * `WG_DNAm_RNA_Upset_Plot.sh`: Use UpSetR to create upset plot showing number of genes genome-wide with differential DNAm and expression across platforms.

* **Supplementary_Data/** (all supplementary files)
  * `Supplementary_Table_1_Lipid_Gene_List.csv`: List of lipid-candidate genes generated using 0_lipid_candidate_gene_extraction.R
  * `Supplementary_Table_2_RNA_DESeq_Lip_Results.xlsx`: Results from DESeq2 differential gene expression analysis in lipid-candidate genes. Combined-sex, sex-stratified, interaction term, and x-chromosomal analyses are reported on each tab.
  * `Supplementary_Table_3_RNA_DESeq_WG_Results.xlsx`: Results from DESeq2 genome-wide differential gene expression analysis . Combined-sex, sex-stratified, interaction term, and x-chromosomal analyses are reported on each tab.
  * `Supplementary_Table_4_Affymetrix_Validation_Lip_Results.xlsx`: Results from differential lipid-candidate gene expression analysis in independent Affymetrix validation cohort. 
  * `Supplementary_Table_5_Affymetrix_Validation_WG_Results.xlsx`: Results from differential gene expression analysis genome-wide in independent Affymetrix validation cohort.
  * `Supplementary_Table_6_Swish_DESeq_rMATS_overlap.xlsx`: Overlap of differentailly expressed genes, transcripts, and alternative splicing identified via Swish, DESeq2, and rMATS.
  * `Supplementary_Table_7_Diff_DNAm_Lip_WG_Results.csv`: Differentially methylated CpGs in lipid-candidate genes and across whole-genome.

* **Results/**
  * `Table_5_Isoform_Switches.csv`: Results of all IsoformSwitchAnalyzeR analyses. File notes name of population group the analysis was run, the name of the upregulated and downregulated transcripts, and any predicted protein consequences.

# Data-Specific Information for: Supplementary_Data/Supplementary_Table_2_RNA_DESeq_Lip_Results.xlsx, Supplementary_Table_3_RNA_DESeq_WG_Results.xlsx
- **Variable List**:
  - `Row.names`: Ensembl ID
  - `baseMean`: Average normalized read counts across all samples
  - `log2FoldChange`: Effect size between PE and Controls
  - `lfcSE`: Standard error for the log2FoldChange
  - `stat`: log2FoldChange/standard error
  - `pvalue`: Unadjusted p-value
  - `padj`: Benjamini-Hochberg FDR adjusted p-value
  - `comparison`: only in interaction term (IT) tab, value denotes which comparison this gene showed significance (`F_PEvsF_Cont`
   = Significant DEG between PE and Controls in Females, `M_PEvsM_Cont` = Significant DEG between PE and Controls in Males, `Interaction_PEvsF_Cont` = Gene has sex-specific differential expression in PE)
  - `gene_symbol`: Gene symbol
  - `seqnames`: Chromosome name
  - `Expression_Status`: Classification of expression change in PE (`Increased_RNA_Expression`
   = log2FC>1,FDR<0.05, `Trending_Towards_Increased_RNA_Expression` = log2FC>0,FDR<0.05, `Decreased_RNA_Expression` = log2FC<-1,FDR<0.05, `Trending_Towards_Decreased_RNA_Expression` = log2FC<0,FDR<0.05, `Not_Biologically_Significant` = FDR>0.05).
  
# Data-Specific Information for: Supplementary_Data/Supplementary_Table_4_Affymetrix_Validation_Lip_Results.xlsx,  Supplementary_Table_5_Affymetrix_Validation_WG_Results.xlsx
- **Variable List**:
  - `probes`:Microarray probe ID
  - `logFC`: Effect size between PE and Controls
  - `AveExpr`: Average log2 expression level across all samples
  - `t`: Moderated t-statistic
  - `P.Value`: Unadjusted p-value
  - `adj.P.Val`: Benjamini-Hochberg FDR adjusted p-value
  - `B`: B-statistic
  - `deltaExprs`: Difference in average expression signal between PE and Controls
  - `SYMBOL`: Gene symbol
  - `GENENAME`: Full gene name
  - `CHR`: Chromosome
  - `ENSEMBL`: Ensembl ID
  - `PMID`:PubMed identifiers
  
# Data-Specific Information for: Supplementary_Data/Supplementary_Table_6_Swish_DESeq_rMATS_overlap

- **Variable List**:
  - `GeneID`: Ensembl ID
  - `X.x`: Artifact from merging columns
  - `Transcript_ID`: Ensembl transcript ID with version number
  - `Tx_ID`: Transcript index from Swish 
  - `Gene_ID.group`: Number corresponding to gene name from Swish
  - `Gene_ID.group_name`: Name corresponding to gene name from Swish (NA if unassigned)
  - `Gene_ID.value`: Full versioned Ensembl gene identifier
  - `Log10_Mean`: log10 mean expression across all samples (Swish)
  - `Stat`: Test statistic for differential transcript expression (Swish)
  - `Log2_FoldChg`: Effect size between PE and Controls (Swish)
  - `P_Value`: Unadjusted p-value (Swish)
  - `Q_Value`: FDR-adjusted q-value (Swish)
  - `Locfdr`: Local FDR estimate for transcript
  - `Gene_Symbol`: Gene symbol (matched from gtf file)
  - `gene_name.x`: Gene symbol (Swish)
  - `X.y`: Artifact from merging columns
  - `baseMean`: Average normalized read counts across all samples (DESeq2)
  - `log2FoldChange`: Effect size between PE and Controls (DESeq2)
  - `lfcSE`: Standard error for the log2FoldChange (DESeq2)
  - `stat`: log2FoldChange/standard error (DESeq2)
  - `pvalue`: Unadjusted p-value (DESeq2)
  - `padj`: Benjamini-Hochberg FDR adjusted p-value (DESeq2)
  - `gene_name.y`: Gene symbol (DESeq2)
  - `seqnames`: Chromosome name (DESeq2)
  - `Expression_Status`: Classification of expression change in PE (`Increased_RNA_Expression`
   = log2FC>1,FDR<0.05, `Trending_Towards_Increased_RNA_Expression` = log2FC>0,FDR<0.05, `Decreased_RNA_Expression` = log2FC<-1,FDR<0.05, `Trending_Towards_Decreased_RNA_Expression` = log2FC<0,FDR<0.05, `Not_Biologically_Significant` = FDR>0.05) (DESeq2)
  - `X`: Artifact of merging columns
  - `geneSymbol`: Gene symbol (rMATS)
  - `chr`: Chromosome (rMATS)
  - `SC_exonStart_0base`/`SC_exonEnd`: Start and end genomic coordinates for Skipped Exon events (rMATS)
  - `SC_FDR`/`SC_IncLevelDifference`: FDR and Inclusion Level Difference for Skipped Exon events (rMATS)
  - `riExonStart_0base`/`riExonEnd`: Genomic coordinates for Retained Intron events (rMATS)
  - `RI_FDR`/`RI_IncLevelDifference`: FDR and Inclusion Level Difference for Retained Intron events (rMATS)
  - `MXE_X1stExonStart_0base`/`MXE_X1stExonEnd`/`MXE_X2ndExonStart_0base`/`MXE_X2ndExonEnd`: Genomic coordinates for Exon 1 and Exon 2 in Mutually Exclusive Exon events (rMATS)
  - `MXE_FDR`/`MXE_IncLevelDifference`: FDR and inclusion level difference for Mutually Exclusive Exon events (rMATS)
  - `A5SS_longExonStart_0base`/`A5SS_longExonEnd`: Genomic coordinates for Alternative 5' Splice Site events (rMATS)
  - `A5SS_FDR`/`A5SS_IncLevelDifference`: FDR and inclusion level difference for Alternative 5' Splice Site events (rMATS)
  - `A3SS_longExonStart_0base`/`A3SS_longExonEnd`: Genomic coordinates for Alternative 3' Splice Site events (rMATS)
  - `A3SS_FDR`/`A3SS_IncLevelDifference`: FDR and inclusion level difference for Alternative 3' Splice Site events (rMATS)


# Data-Specific Information for: Supplementary_Data/Supplementary_Table_7_Diff_DNAm_Lip_WG_Results

- **Variable List**:
  - `X`: Artifact of merging columns
  - `probes`: Probe ID
  - `logFC`: Log2foldchange of methylation levels between PE and controls
  - `AveExpr`: Average signal intensity across all samples
  - `t`: Moderated t-statistic
  - `P.Value`: Unadjusted p-value
  - `adj.P.Val`: Benjamini-Hochberg FDR adjusted p-value
  - `B`: B-statistic
  - `deltaB`: Change in beta value (methylation level difference between groups)
  - `gene`: Gene symbol
  - `chr`: Chromosome name
  - `position`: Genomic coordinate position of the probe
  - `region_overlap`: Genomic region or feature the probe overlaps with (e.g., promoter, gene body, enhancer)
  - `Closest_TSS_gene_name`: Gene name of the closest Transcription Start Site (TSS)
  - `Closest_TSS_Transcript`: Transcript ID of the closest TSS
  - `Closest_TSS_Pos`: Genomic position of the closest TSS
  - `Distance_Closest_TSS`: Genomic distance from the probe to the closest TSS
  - `gene_id`: Ensembl ID
  - `overlap`: Overlap between gene name and closest TSS
  - `diffmethylation`: Classification of differential methylation in PE (`Increased_RNA_Expression`
   = deltaB>0.05,FDR<0.05, `Trending_Towards_Increased_RNA_Expression` = deltaB>0,FDR<0.05, `Decreased_RNA_Expression` = deltaB<-0.05,FDR<0.05, `Trending_Towards_Decreased_RNA_Expression` = deltaB<0,FDR<0.05, `Not_Biologically_Significant` = FDR>0.05).

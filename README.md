This readme file was generated on 2025-02-18 by Kriesha Eyer.

# General Information

Title of Dataset: 2024_PE_lipid_DNA_methylation

Author/Principal Investigator Information
Name: Samantha Wilson
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

Date of data collection: 2024-10-09 - 2025-02-28 

Geographic location of data collection: Hamilton, Ontario, Canada

# Sharing/Access Information

Was data derived from another source?

### Lipid Candidate Gene Search: 
- GAD Gene-Disease Associations Dataset - https://maayanlab.cloud/Harmonizome/dataset/GAD+Gene-Disease+Associations
- GO Biological Process Annotations 2023 Dataset - https://maayanlab.cloud/Harmonizome/dataset/GO+Biological+Process+Annotations+2023
-  GO Cellular Component Annotations 2023 Dataset - https://maayanlab.cloud/Harmonizome/dataset/GO+Cellular+Component+Annotations+2023
- Reactome Pathways 2024 Dataset - https://maayanlab.cloud/Harmonizome/dataset/Reactome+Pathways+2024
  
### Gene Expression Omnibus: 
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125605
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100197
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE98224

# Data & File Overview

All project scripts and plots for Differential DNA Methylation of Lipid Candidate Genes in Preeclampsia. Order of read me is the order scipts should be run.

### lipid_candidate_gene_extraction.R

Identify lipid-related genes from Harmonizome data sets based on key lipid-related terms.

### meta_data_extraction.R

Extract clinical and technical data for all samples in GEO data sets using GEOquery.

### fetal_sex_prediction.R

Use minfi package to estimate fetal sex for all GEO data sets as quality control measure and to add fetal sex as a variable. 

### normalization_probefiltering.R

1. Load data sets into R and exclude samples that were not control or preeclamptic patients in studies (ex. IUGR samples).
2. Normalize data using adjustedFunnorm.
3. Exclude probes with bad detection P-values, missing beta-values, SNP probes, cross-hybridizing probes and non-variable placental probes. 

### methylation_analysis.R

1. Stratifying groups by sex for methylation analysis.
2. Use delta beta values between control and preeclamptic patients to identify significant differentially methylated lipid-related genes.
3. Use M-values to perform linear modeling for visualization of results using volcano plots.

### density_plots_beta_values.R

Plot the adjusted P-values, P-values, delta-beta-values and beta-values to identify whether there are any outlier samples. 

### pathway_enrichment_analysis.R

Using ErmineJ, identify the most common differentially methylated pathways in the data set.




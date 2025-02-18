# 2024_PE_lipid_DNA_methylation

All project scripts and plots for Lipid Gene DNA methyation Analysis in Preeclamptic Patients - order of read me is the order scipts should be run

# lipid_candidate_gene_extraction.R

Identify lipid-related genes from Harmonizome data sets based on key lipid-related terms.

# meta_data_extraction.R

Extract clinical and technical data for all samples in GEO data sets using GEOquery.

# fetal_sex_prediction.R

Use minfi package to estimate fetal sex for all GEO data sets as quality control measure and to add fetal sex as a variable. 

# normalization_probefiltering.R

1. Load data sets into R and exclude samples that were not control or preeclamptic patients in studies (ex. IUGR samples).
2. Normalize data using adjustedFunnorm.
3. Exclude probes with bad detection P-values, missing beta-values, SNP probes, cross-hybridizing probes and non-variable placental probes. 

# methylation_analysis.R

1. Stratifying groups by sex for methylation analysis.
2. Use delta beta values between control and preeclamptic patients to identify significant differentially methylated lipid-related genes.
3. Use M-values to perform linear modeling for visualization of results using volcano plots.

# density_plots_beta_values.R

Plot the adjusted P-values, P-values, delta-beta-values and beta-values to identify whether there are any outlier samples. 

# pathway_enrichment_analysis.R

Using ErmineJ, identify the most common differentially methylated pathways in the data set.




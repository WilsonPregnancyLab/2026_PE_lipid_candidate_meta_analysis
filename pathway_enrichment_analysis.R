# Part 1: NCBI Reference Subsetting
# the goal here is to add the NCBI annotations to the lipid genes (First, I'm going to do it for all of the lipid genes I tested - in whole_pop, M and F, then only for the significant ones in each group)

 if(!file.exists('Generic_human_ncbiIds_noParents.an.txt.gz')){
    system('wget https://gemma.msl.ubc.ca/annots/Generic_human_ncbiIds_noParents.an.txt.gz')}
NCBI <- read.table('Generic_human_ncbiIds_noParents.an.txt.gz', sep = '\t', header = T, quote = "")

# Whole Population Autosomes
placmet_wholepop_auto <- read.csv("placmet_wholepop_auto.csv")
placmet_wholepop_auto_pval <- placmet_wholepop_auto[, c("Closest_TSS_gene_name", "adj.P.Val")] # 46613 probes
colnames(placmet_wholepop_auto_pval)[colnames(placmet_wholepop_auto_pval) == "Closest_TSS_gene_name"] <- "GeneSymbols"
placmet_wholepop_auto_reference <- NCBI[NCBI$GeneSymbols %in% placmet_wholepop_auto_pval$GeneSymbols, ] #2845 genes
placmet_wholepop_auto_pval <- placmet_wholepop_auto_pval[!duplicated(placmet_wholepop_auto_pval$GeneSymbols),] #2887 lipid genes investigated
write.csv(placmet_wholepop_auto_pval, "placmet_wholepop_auto_pval.csv")
write.csv(placmet_wholepop_auto_reference, "placmet_wholepop_auto_NCBI_reference.csv")

# Female Autosomes
placmet_F_auto <- read.csv("placmet_F_fulldata_auto.csv")
placmet_F_auto_pval <- placmet_F_auto[, c("Closest_TSS_gene_name", "adj.P.Val")] # 46613 probes
colnames(placmet_F_auto_pval)[colnames(placmet_F_auto_pval) == "Closest_TSS_gene_name"] <- "GeneSymbols"
write.csv(placmet_F_auto_pval, "placmet_F_auto_pval.csv")

# Male Autosomes
placmet_M_auto <- read.csv("placmet_M_fulldata_auto.csv")
placmet_M_auto_pval <- placmet_M_auto[, c("Closest_TSS_gene_name", "adj.P.Val")] # 46613 probes
colnames(placmet_M_auto_pval)[colnames(placmet_M_auto_pval) == "Closest_TSS_gene_name"] <- "GeneSymbols"
write.csv(placmet_M_auto_pval, "placmet_M_auto_pval.csv")


# Part 2: Pathway Enrichment Analysis
devtools::install_github('PavlidisLab/ermineR')
library(ermineR)
library(dplyr)
library(ggplot2)
install.packages("rJava")

Sys.setenv('JAVA_HOME' = '/usr/lib/jvm/java-21-openjdk-amd64/')

placmet_wholepop_auto_pval <- read.csv("placmet_wholepop_auto_pval.csv")
placmet_wholepop_auto_reference <- read.csv("placmet_wholepop_auto_NCBI_reference.csv")
placmet_wholepop_F_pval <- read.csv("placmet_F_auto_pval.csv")
placmet_wholepop_M_pval <- read.csv("placmet_M_auto_pval.csv")

all_wholepop_auto_pr_out <- precRecall(annotation = placmet_wholepop_auto_reference, 
                    scores = placmet_wholepop_auto_pval,
                    scoreColumn = 2,
                    logTrans = T,
                    bigIsBetter = F,
                    iterations = 10000)

all_wholepop_auto_pathway_analysis <- as.data.frame(all_wholepop_auto_pr_out$results)

sig_pathways_wholepop_auto <- all_wholepop_auto_pathway_analysis[all_wholepop_auto_pathway_analysis$CorrectedPvalue < 0.05 & all_wholepop_auto_pathway_analysis$Multifunctionality < 0.5,] #2
sig_pathways_wholepop_auto <- sig_pathways_wholepop_auto %>% distinct(GeneMembers, .keep_all = T)
write.csv(sig_pathways_wholepop_auto, '2025_wholepop_auto_enriched_lipid_pathways.csv')
# 2887 lipid genes & filtered reference: there were 656 pathways in the analysis, 2 sig enriched pathways

all_F_auto_pr_out <- precRecall(annotation = placmet_wholepop_auto_reference, 
                    scores = placmet_F_auto_pval,
                    scoreColumn = 2,
                    logTrans = T,
                    bigIsBetter = F,
                    iterations = 10000)

all_F_auto_pathway_analysis <- as.data.frame(all_F_auto_pr_out$results)

sig_pathways_F_auto <- all_F_auto_pathway_analysis[all_F_auto_pathway_analysis$CorrectedPvalue < 0.05 & all_F_auto_pathway_analysis$Multifunctionality < 0.5,] #2
sig_pathways_F_auto <- sig_pathways_F_auto %>% distinct(GeneMembers, .keep_all = T)
write.csv(sig_pathways_F_auto, '2025_F_auto_enriched_lipid_pathways.csv')
# 2887 lipid genes & filtered reference: there were 656 pathways in the analysis, 0 sig enriched pathways

all_M_auto_pr_out <- precRecall(annotation = placmet_wholepop_auto_reference, 
                    scores = placmet_M_auto_pval,
                    scoreColumn = 2,
                    logTrans = T,
                    bigIsBetter = F,
                    iterations = 10000)

all_M_auto_pathway_analysis <- as.data.frame(all_M_auto_pr_out$results)

sig_pathways_M_auto <- all_M_auto_pathway_analysis[all_M_auto_pathway_analysis$CorrectedPvalue < 0.05 & all_M_auto_pathway_analysis$Multifunctionality < 0.5,] #2
sig_pathways_M_auto <- sig_pathways_M_auto %>% distinct(GeneMembers, .keep_all = T)
write.csv(sig_pathways_M_auto, '2025_M_auto_enriched_lipid_pathways.csv')
# 2887 lipid genes & filtered reference: there were 1764 pathways in the analysis, 0 sig enriched pathways

# Plotting

ord_sig_pathways_wholepop_auto <- sig_pathways_wholepop_auto [order(sig_pathways_wholepop_auto$CorrectedPvalue, decreasing = F),]
rownames(ord_sig_pathways_wholepop_auto) <- NULL

ord_sig_pathways_wholepop_auto$Name <- factor(ord_sig_pathways_wholepop_auto$Name, levels = rev(ord_sig_pathways_wholepop_auto$Name[order(ord_sig_pathways_wholepop_auto$CorrectedPvalue)]))

ord_sig_pathways_wholepop_auto$Label <- paste0(ord_sig_pathways_wholepop_auto$Name, '(', ord_sig_pathways_wholepop_auto$ID, ')')
ord_sig_pathways_wholepop_auto$Label <- factor(ord_sig_pathways_wholepop_auto$Label, levels = rev(ord_sig_pathways_wholepop_auto$Label[order(ord_sig_pathways_wholepop_auto$CorrectedPvalue)]))

pathway_wholepop_auto_plot <- ggplot(ord_sig_pathways_wholepop_auto[1:2,], aes(x = Label, y = -log(CorrectedPvalue), fill = Multifunctionality)) +
    geom_bar(stat = 'identity') +
    # geom_hline(yintercept = -log(0.05), col = '#a30000', linetype = 'dashed') +
    xlab('') +
    ylab('-log(FDR)') +
    scale_fill_gradient2(
        low = "#7727e0",
        mid = "#ffa837",
        high = "#ffa837", 
        midpoint = 0.5,
        limits = c(0, 0.5),
        oob = scales::squish
    ) +
    coord_flip() +
    theme_light() +
    theme(legend.position = 'bottom',
    legend.justification = c(3.5,0),
    text = element_text(size = 15),
    axis.text.x = element_text(color = 'black'),
    axis.text.y = element_text(color = 'black'),
    legend.spacing.y = unit(0.3, 'cm'),
    legend.title = element_text(margin = margin(r = 40, b = 10))) +
    guides(fill = guide_colorbar(barwidth = 5, barheight = 0.7))

png(filename = './2025_Whole_Pop_Auto_Lipid_DMR_Pathways.png', width = 15, height = 10, units = 'in', res = 300)
pathway_wholepop_auto_plot
dev.off()



wholepop_auto_enriched_lipid_pathways















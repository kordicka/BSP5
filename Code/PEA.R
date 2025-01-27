# Install necessary libraries if not already installed
# install.packages("BiocManager")
# BiocManager::install("clusterProfiler")
# BiocManager::install("org.Rn.eg.db")
# BiocManager::install("enrichplot")
# BiocManager::install("DOSE")
# install.packages("tidyverse")
# install.packages("pheatmap")
# install.packages("plotly")


# Load required libraries
library(clusterProfiler)
library(org.Rn.eg.db)
library(enrichplot)
library(ggplot2)
library(DOSE)
library(plotly)

# Read in the differential expression results
diff_expr <- read.table("/home/nora/Semester5/BSP5/differential_expression_results.tsv", header=TRUE, sep="\t")

# Extract significant genes (adjusted p-value < 0.05)
sig_genes <- subset(diff_expr, adj.P.Val < 0.05)$Gene.symbol

# Convert gene symbols to Entrez IDs
gene_ids <- bitr(sig_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Rn.eg.db)

# Perform GO enrichment analysis
go_enrichBP <- enrichGO(gene          = gene_ids$ENTREZID,
                      OrgDb         = org.Rn.eg.db,
                      ont           = "BP", 
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05)

go_enrichCC <- enrichGO(gene          = gene_ids$ENTREZID,
                      OrgDb         = org.Rn.eg.db,
                      ont           = "CC", 
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05)

go_enrichMF <- enrichGO(gene          = gene_ids$ENTREZID,
                      OrgDb         = org.Rn.eg.db,
                      ont           = "MF", 
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05)

# Perform KEGG pathway enrichment analysis
kegg_enrich <- enrichKEGG(gene         = gene_ids$ENTREZID,
                          organism     = "rno", 
                          pvalueCutoff = 0.05)

# Ensure the pairwise similarity matrix for tree and enrichment map plots
go_enrichBP <- pairwise_termsim(go_enrichBP)
kegg_enrich <- pairwise_termsim(kegg_enrich)

#################PLOTS
# Interactive Dot Plot for GO
go_dotBP <- dotplot(go_enrichBP, showCategory=10) + ggtitle("GO Enrichment (BP) - Dot Plot")
ggsave("go_dotBP.png", plot = go_dotBP, width = 8, height = 6, dpi = 300)
ggplotly(go_dotBP)

# Interactive Dot Plot for KEGG
kegg_dot <- dotplot(kegg_enrich, showCategory=10) + ggtitle("KEGG Enrichment - Dot Plot")
ggsave("kegg_dot.png", plot = kegg_dot, width = 8, height = 6, dpi = 300)
ggplotly(kegg_dot)

# Interactive Heatmap for GO
go_heatBP <- heatplot(go_enrichBP, showCategory=10) + ggtitle("GO Enrichment(BP) - Heatmap")
ggplotly(go_heatBP)

# Interactive Heatmap for KEGG
kegg_heat <- heatplot(kegg_enrich, showCategory=10) + ggtitle("KEGG Enrichment - Heatmap")
ggplotly(kegg_heat)

# Static Network-Based Plots
# Set PDF device for saving static plots
pdf("/home/nora/Semester5/BSP5/static_network_plots.pdf", width = 8, height = 6)

# Cnet Plot (Network Plot) for GO
cnetplot(go_enrichBP, showCategory=5) + ggtitle("GO Enrichment (BP) Network Plot")

# Cnet Plot for KEGG
cnetplot(kegg_enrich, showCategory=5) + ggtitle("KEGG Enrichment Network Plot")

# Tree Plot for GO
treeplot(go_enrichBP) + ggtitle("GO Enrichment (BP) Tree Plot")

# Tree Plot for KEGG
treeplot(kegg_enrich) + ggtitle("KEGG Enrichment Tree Plot")

# Enrichment Map for GO
emapplot(go_enrichBP, showCategory = 10) + ggtitle("GO Enrichment (BP) Map")

# Enrichment Map for KEGG
emapplot(kegg_enrich, showCategory = 10) + ggtitle("KEGG Enrichment Map")

###############SORTING
# Sort pathways by adjusted p-value
go_results <- go_enrichBP@result
sorted_go_results <- go_results[order(go_results$p.adjust), ]

# Filter for pathways with p.adjust < 0.05
significant_go_results <- subset(sorted_go_results, p.adjust < 0.05)

# Create the core_genes column
significant_go_results$core_genes <- sapply(significant_go_results$geneID, function(g) paste(strsplit(g, "/")[[1]], collapse=", "))

# Display the top 10 most significant pathways
head(significant_go_results[, c("Description", "p.adjust", "core_genes")], 10)

# Sort KEGG pathways by adjusted p-value
kegg_results <- kegg_enrich@result
sorted_kegg_results <- kegg_results[order(kegg_results$p.adjust), ]

# Filter for pathways with p.adjust < 0.05
significant_kegg_results <- subset(sorted_kegg_results, p.adjust < 0.05)

# Create the core_genes column
significant_kegg_results$core_genes <- sapply(significant_kegg_results$geneID, function(g) paste(strsplit(g, "/")[[1]], collapse=", "))

# Display the top 10 most significant KEGG pathways
head(significant_kegg_results[, c("Description", "p.adjust", "core_genes")], 10)

# Close PDF device
dev.off()

##########CREATON OF FILES FOR NETWORK ANALYSIS
# 1. File with gene symbols
gene_symbols <- diff_expr["Gene.symbol"]
write.table(gene_symbols, "/home/nora/Semester5/BSP5/gene_symbols.tsv", sep="\t", row.names=FALSE, quote=FALSE)

# 2. File with gene symbols and log fold changes
gene_symbols_logfc <- diff_expr[, c("Gene.symbol", "logFC")]
write.table(gene_symbols_logfc, "/home/nora/Semester5/BSP5/gene_symbols_logfc.tsv", sep="\t", row.names=FALSE, quote=FALSE)

# 3. File with gene symbols and negative log fold changes
gene_symbols_neg_logfc <- gene_symbols_logfc
gene_symbols_neg_logfc$logFC <- -gene_symbols_neg_logfc$logFC
write.table(gene_symbols_neg_logfc, "/home/nora/Semester5/BSP5/gene_symbols_neg_logfc.tsv", sep="\t", row.names=FALSE, quote=FALSE)

# 4. File with gene symbols and 0/1 for regulation (downregulated = 0, upregulated = 1)
gene_symbols_up_down <- gene_symbols_logfc
gene_symbols_up_down$Regulation <- ifelse(gene_symbols_up_down$logFC > 0, 1, 0)
gene_symbols_up_down <- gene_symbols_up_down[, c("Gene.symbol", "Regulation")]
write.table(gene_symbols_up_down, "/home/nora/Semester5/BSP5/gene_symbols_up_down.tsv", sep="\t", row.names=FALSE, quote=FALSE)

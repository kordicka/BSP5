# Load necessary libraries
library(GEOquery)
library(limma)
library(umap)
library(ggrepel)
library(ggplot2)

# Load and prepare GEO data
gset <- getGEO("GSE139438", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]
fvarLabels(gset) <- make.names(fvarLabels(gset))

# Group assignment
gsms <- "111111000000"
gs <- factor(strsplit(gsms, split = "")[[1]], labels = c("treated", "untreated"))

# Log2 transformation if required
ex <- exprs(gset)
if (any(quantile(ex, 0.99) > 100)) { 
  ex[ex <= 0] <- NaN
  exprs(gset) <- log2(ex)
}

# Set up design matrix and fit model
design <- model.matrix(~ gs + 0)
colnames(design) <- levels(gs)
fit <- eBayes(contrasts.fit(lmFit(gset, design), makeContrasts("treated-untreated", levels = design)))

# Extract and save top significant genes
tT <- topTable(fit, adjust = "fdr", sort.by = "B", number = 250)
write.table(tT[, c("ID", "adj.P.Val", "P.Value", "t", "B", "logFC", "Gene.symbol", "Gene.title")], 
            file = stdout(), row.names = FALSE, sep = "\t")

# Visualization and QC
hist(tT$adj.P.Val, col = "grey", main = "P-adj value distribution", xlab = "P-adj")
dT <- decideTests(fit, adjust.method = "fdr", p.value = 0.05, lfc = 0)
vennDiagram(dT, circle.col = palette())
qqt(fit$t[!is.na(fit$F)], fit$df.total[!is.na(fit$F)], main = "Moderated t statistic")
volcanoplot(fit, coef = 1, main = "Volcano plot", pch = 20, highlight = sum(dT != 0))
plotMD(fit, column = 1, status = dT[, 1], legend = FALSE, pch = 20)

write.table(tT, file="differential_expression_results.tsv", row.names=FALSE, sep="\t", quote=FALSE)

if (!require("BiocManager")) 
  install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Sc.sgd.db", "enrichplot"))
library(DESeq2)
library(tximport)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(VennDiagram)
library(clusterProfiler)
library(org.Sc.sgd.db)
library(enrichplot)

#Setting the working directory
#Loading all the RSEM files needed for analysis
dir <- "~/ubunututo windows/rsem_output"
samples <- read.csv("~/ubunututo windows/metadata.csv")
files <- file.path(dir, paste0(samples$sample, ".genes.results"))
names(files) <- samples$sample
#Import data using trixmport
txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
#Checking for the top genes
count_matrix <- txi.rsem$counts
head(count_matrix)
#Creating the DESeq2 object
dds <- DESeqDataSetFromTximport(txi.rsem, colData = samples, design = ~ stage)
# The below function performs dispersion estimate
dds <- DESeq(dds)
colnames(colData(dds))
levels(dds$stage) #the three stages of the biofilm
resultsNames(dds) #verify how the model is structured

#1. Relationship between samples using the developmental stages of Biofilm
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = "stage") +
  ggtitle("PCA of Biofilm Development Stages")
#Extract differential expression results for two structures using the model coefficient
#Mature vs Early Biofilm
res_mature_vs_early <- results(dds, name="stage_.mature_vs_.early", alpha=0.05)
#Thin vs Early Biofilm
res_thin_vs_early <- results(dds, name="stage_.thin_vs_.early", alpha=0.05)

# Look at the statistical overview of the comparison
summary(res_mature_vs_early) 
summary(res_thin_vs_early)

#genes most upregulated in Mature vs Early using log2FoldChange
res_up <- res_mature_vs_early[order(res_mature_vs_early$log2FoldChange, decreasing=TRUE), ]
head(res_up, 20)
#genes most downregulated
res_down <- res_mature_vs_early[order(res_mature_vs_early$log2FoldChange), ]
head(res_down, 20)
#ordering genes by padj(statistically significant)
res_mature_vs_early <- res_mature_vs_early[order(res_mature_vs_early$padj), ]
res_thin_vs_early <- res_thin_vs_early[order(res_thin_vs_early$padj), ]
head(res_thin_vs_early)
head(res_mature_vs_early)

#Total count of significant genes
sum(res_mature_vs_early$padj < 0.05, na.rm=TRUE)
sum(res_thin_vs_early$padj < 0.05, na.rm=TRUE)

#MA plots for the two comparison
#For the MA plots, created my object as a data frame, introduced padj<0.05
#2a
res_dfME <- as.data.frame(res_mature_vs_early)
res_dfME$Significance <- ifelse(res_dfME$padj < 0.05, 
                                 "Significant", 
                                 "Not Significant")
ggplot(res_dfME, aes(x = baseMean, y = log2FoldChange, color = Significance)) +
  geom_point(alpha = 0.6, size = 1.0) +
  scale_x_log10() +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "blue")) +
  geom_hline(yintercept = 0, linetype = "solid") + ylim(-5, 5) +
  labs(title = "MA Plot: Mature vs Early", x = "Mean of Normalized Counts (log scale)",
       y = "Log2 Fold Change (Mature vs Early)") +
  theme_minimal(base_size = 10) + theme(legend.title = element_blank())
#b
res_dfTE <- as.data.frame(res_thin_vs_early)
res_dfTE$Significance <- ifelse(res_df_TE$padj < 0.05, 
                                 "Significant", 
                                 "Not Significant")
ggplot(res_dfTE, aes(x = baseMean, y = log2FoldChange, color = Significance)) +
  geom_point(alpha = 0.6, size = 1.0) +
  scale_x_log10() +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "blue")) +
  geom_hline(yintercept = 0, linetype = "solid") + ylim(-5, 5) +
  labs(title = "MA Plot: Thin vs Early", x = "Mean of Normalized Counts (log scale)",
       y = "Log2 Fold Change (Thin vs Early)") +
  theme_minimal(base_size = 10) + theme(legend.title = element_blank())
#3 Volcano plot using mature vs early comparison
res_df <- as.data.frame(res_mature_vs_early)
res_df$Significance <- "NS"#Added a new column so as to distinguish the up, down regulated genes
res_df$Significance[res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "Up"
res_df$Significance[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Down"
colnames(res_df) 
# plot
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
  geom_point(alpha=0.4, size=1.2) +
  scale_color_manual(values=c("Up"="red", "Down"="blue", "NS"="grey60")) +
  theme_minimal() +
  labs(title = "DESeq2 Volcano Plot-ME", x = "log2 Fold Change", y = "-log10 Adjusted p-value")

#4 Get names of top 10 genes by adjusted p-value
top_genes <- head(order(res_mature_vs_early$padj), 10)
mat <- assay(vsd)[top_genes, ]
mat <- mat - rowMeans(mat) # Center the data
# Visualizing the genes
pheatmap(mat, annotation_col=as.data.frame(colData(dds)["stage"]),
         main="Top 10 Differentially Expressed Genes")

#5
# Detecting and filtering the significant genes in mature vs early biofilm and thin vs early
#based on magnitude and statistical significance. Get significant gene names for both comparisons
# Using padj < 0.05  and log2fold change as the threshold
sig_genesME <- rownames(subset(res_mature_vs_early, padj < 0.05 & abs(log2FoldChange) > 1))
sig_genesTE <- rownames(subset(res_thin_vs_early, padj < 0.05 & abs(log2FoldChange) > 1))
length(intersect(sig_genesME, sig_genesTE))
#Generate the Venn Diagram
venn.plot <- venn.diagram(
  x = list(Thin_vs_Early = sig_genesTE, Mature_vs_Early = sig_genesME),
  filename = NULL, 
  fill = c("blue", "green"),
  alpha = 0.6, cex = 1.0, cat.cex = 1.0, main = "Overlap of Differentially Expressed Genes")
grid.draw(venn.plot)

#6. now let us Perform GO Enrichment for each comparisons

# Combine the gene lists into a named list
genescombined <- list(ThinEarly= sig_genesTE, MatureEarly = sig_genesME)

# Run the comparison across Biological Process (BP)
compGO <- compareCluster(geneCluster = genescombined,
                         fun = "enrichGO", 
                         OrgDb = org.Sc.sgd.db, 
                         keyType = "ENSEMBL",
                         ont = "BP", 
                         pAdjustMethod = "BH")
head(compGO)
# Visualize the comparison
dotplot(compGO, showCategory = 3) + 
  ggtitle("Biofilm stages")

#end
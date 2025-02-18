#### Libraries DONE ####
library(DESeq2)
library(dplyr)
library(NMF)
library(ggplot2)
library(clusterProfiler)

#### Import Data ####
load("D:/Assignments/CB1/RNA_Assignment/gene.Rdata")
load("D:/Assignments/CB1/RNA_Assignment/norm.RData")

# Convert Entrez ID to Gene Symbols
for (i in 1:nrow(Counts)) {
  # Extract the current GeneID from Counts
  current_gene_id <- rownames(Counts)[i]
  
  gene_symbol <- GeneInfo[GeneInfo$GeneID == as.integer(rownames(Counts)[i]),2]
  
  # Assign the GeneSymbol as the new row name if it exists
  if (length(gene_symbol) >0){
  rownames(Counts)[i] <- gene_symbol
  }
  print(i)
}

# Create names of the samples
sampleNames <- c()
nameAPR <- c()
nameDMSO <- c()

for (i in 1:3){
  nameAPR <- paste ("APR_", i, sep = "")
  sampleNames <- c(sampleNames, nameAPR)
  nameDMSO <- paste ("DMSO_", i, sep = "")
  sampleNames <- c(sampleNames, nameDMSO)
}

# Create the metadata for running the DGE
colnames(Counts) <- sort(sampleNames)
Counts <- Counts[,c(1,2,3,4,6)]
metadata <- data.frame(
  condition = factor(c("Treated", "Treated", "Treated", "Control", "Control")),
  row.names = colnames(Counts)
)

#### Process Counts Data ####
# Create DESeq object
deseqRaw <- DESeqDataSetFromMatrix(countData = Counts, 
                                   colData = metadata,
                                   design = ~ c("Treated", "Treated", "Treated", "Control", "Control"))
# Remove reads with no expression and normalise
deseqRaw <- deseqRaw[ rowSums(counts(deseqRaw)) > 0, ]
deseqNorm <- estimateSizeFactors(deseqRaw)

# Transform the data to look for normality and outliers
countsNorm <- counts(deseqNorm, normalized = T)
deseqRlog <- rlog(deseqNorm, blind = T)
normalizedCountsRlog <- assay(deseqRlog)
logNormCounts <- log2(countsNorm + 1)

# Hierarchical clustering object
distanceRlog <- as.dist(1- cor(normalizedCountsRlog, method = "spearman"))

# Make the Control condition the baseline
colData(deseqRaw)$condition <- relevel(colData(deseqRaw)$condition, "Control")

# Create DESeq matrix and store results
deseqRawMatrix <- DESeq(deseqRaw)
deseqDgeResults <- results(deseqRawMatrix, independentFiltering = T, alpha = 0.05)

#### Data Summaries #####

# Summaries data from the raw matrix
colData (deseqRaw) %>% head
assay (deseqRaw, "counts") %>% head
rowData (deseqRaw) %>% head
counts (deseqRaw) %>% str

# These values should all be true
colSums (counts(deseqRaw)) == colSums (Counts)

sizeFactors(deseqNorm)
colData (deseqNorm) %>% head 

#### Visualise Data #####
# Visualise the normalisation and trnasformation
plot(logNormCounts[,1:2], cex = 1, main = "Normalized and Log2 Transformed")
plot(normalizedCountsRlog[,1:2], cex = 1, main = "Regularised Log-Transformed (Rlog)")

# Create a hierachical clustering plot and PCA
plot(hclust(distanceRlog), labels = colnames(normalizedCountsRlog))
pcaPlot <- plotPCA(deseqRlog)
pcaPlot

# Create histograms and MA plots
hist(deseqDgeResults$padj, main = "Histogram of padj from the DGE results")
DESeq2::plotMA(deseqDgeResults[deseqDgeResults$baseMean > 15,], 
               alpha = 0.05, ylim = c(-3,5), main = "MA-Plot of DEgenes")

# Summaries results
summary(deseqDgeResults)
head (deseqDgeResults)
table (deseqDgeResults$padj < 0.05)

# Create heatmap variables and heatmap
heatmapGenes <- logNormCounts[!is.na(deseqDgeResults$padj) 
                              & deseqDgeResults$baseMean > 15
                              & deseqDgeResults$padj < 0.05, ]

orderedGenes <- deseqDgeResults[rownames(heatmapGenes),]
orderedGenes <- orderedGenes[order(orderedGenes$log2FoldChange, decreasing = T),]
heatmapGenes <- heatmapGenes[rownames(orderedGenes),]
top20Genes <- heatmapGenes[c(1:10, (nrow(heatmapGenes)-9):nrow(heatmapGenes)), ]
aheatmap(heatmapGenes, Rowv = NA, Colv = NA)
aheatmap(top20Genes, Rowv = NA, Colv = NA)
rownames(top20Genes)

# Create a volcano plot which labels significant and non significant genes different colours
volcano <- deseqDgeResults[!is.na(deseqDgeResults$padj),]
volcano$significance <- ifelse(volcano$padj < 0.05, "Significant", "Not Significant")

ggplot(volcano, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("blue", "red")) +
  theme_minimal() +
  labs(
    title = "Volcano Plot",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value") +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank())

#### Explore DGE ####
# This section manipulates the DEGs to view most significant results
sortedDeseqDgeResults <- deseqDgeResults[order(deseqDgeResults$padj), ]
degenesPadjDeseq <- subset(sortedDeseqDgeResults, padj < 0.05)
degenesPadjLogDeseq <- subset(sortedDeseqDgeResults,
                              padj < 0.05 &
                                abs(log2FoldChange) >= 1)

sortLog <- degenesPadjLogDeseq[order(degenesPadjLogDeseq$log2FoldChange),]
overTop5 <- tail(sortLog, 6)
underTop5 <- head(sortLog, 6)

rownames(overTop5)
rownames(underTop5)

# Visualise the change in expression for a particular gene, in this case TP53
plotCounts(dds = deseqRaw,
           gene = "TP53",
           normalized = T, transform = F)

deseqDgeResults["TP53",]

#### Geneset Enrichment #####
# Create table for GSEA
testGSEA <- deseqDgeResults[order(-deseqDgeResults$stat),]
testGSEA <- data.frame(rownames(testGSEA), testGSEA[,"stat"])

write.table(testGSEA, "ranking.rnk",
            sep = "\t", 
            col.names = F,
            row.names = F)

#### Pathway Enrichment #####
design <- model.matrix(~ 0+factor(c(rep(1,3),rep(2,3))))
colnames(design) <- c("DMSO", "APR")
contrast.matrix <- makeContrasts(APR-DMSO, levels=design)
v <- voom(y, design)
fit <- lmFit(v, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit3 <- eBayes(fit2)
topTable(fit3)
topTable(fit3,number=1000,lfc=1.2,p.value=0.05)

de <- topTreat(fit3, number=1000,p.value=0.05)[,1]
yy <- enrichKEGG(de, pvalueCutoff=0.05)
head(yy)

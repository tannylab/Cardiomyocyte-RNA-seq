library(readr)
library(tximport)
library(DESeq2)
library(IHW)
library(pheatmap)

# List StringTie output files for all samples
# All files should be in same directory
files <- list.files(here::here(''), pattern = "t_data.ctab", recursive = TRUE, full.names = TRUE)

# Read in StringTie files with tximport
tmp <- read_tsv(files[1])
tx2gene <- tmp[, c("t_name", "gene_name")]
txi <- tximport(files, type = "stringtie", tx2gene = tx2gene)

# Extract sample information from file paths
sampleID <- unlist(lapply(strsplit(files, "/"), "[", 11))
agonist.antagonist <- lapply(strsplit(sampleID, "_"), "[", c(2,3))
agonist.antagonist <- lapply(agonist.antagonist, function(x) paste(x[1], x[2], sep = '.'))
agonist.antagonist <- unlist(agonist.antagonist)
agonist.antagonist <- gsub("-", "", agonist.antagonist)
batch <- as.factor(unlist(lapply(strsplit(sampleID, "_"), "[", 4)))

# Put sample information into dataframe for DESeq2
sample.table <- data.frame(sampleID=sampleID,
                           agonist.antagonist=agonist.antagonist,
                           batch=batch)

# Relevel agonist.antagonist 
sample.table$agonist.antagonist <- relevel(factor(sample.table$agonist.antagonist), ref="Acetic.DMSO")

# Create DESeq Object including batch in design
dds <- DESeqDataSetFromTximport(txi,
                                sample.table,
                                ~ agonist.antagonist + batch)

# Remove low abundance rows
dds <- dds[ rowSums(counts(dds)) > 1, ]

# Run analysis
dds <- DESeq(dds)

# Obtain log2 Fold Change for required comparisons
resIHW.PEonly <- results(dds, filterFun=ihw, alpha=0.05, contrast=c("agonist.antagonist", "PE.DMSO", "Acetic.DMSO"))
resIHW.pejq1.pedmso <- results(dds, filterFun=ihw, alpha=0.05, contrast=c("agonist.antagonist", "PE.JQ1", "PE.DMSO"))

resIHW.ET1only <- results(dds, filterFun=ihw, alpha=0.05, contrast=c("agonist.antagonist", "ET1.DMSO", "Acetic.DMSO"))
resIHW.et1jq1.et1dmso <- results(dds, filterFun=ihw, alpha=0.05, contrast=c("agonist.antagonist", "ET1.JQ1", "ET1.DMSO"))

differential_genes <- function(agonist.only, jq1.agonist, agonist.name) {
  # Extract genes upregulated by PE
  Up <- agonist.only[(agonist.only$log2FoldChange > 1) & (agonist.only$padj < 0.05),]
  Up <- unlist(rownames(Up))

  # Extract genes downregulated by PE
  Down <- agonist.only[(agonist.only$log2FoldChange < -1) & (agonist.only$padj < 0.05),]
  Down <- unlist(rownames(Down))

  # Extract PE upregulated genes attenuated by JQ1
  Up_JQ1 <- jq1.agonist
  Up_JQ1$Gene <- rownames(jq1.agonist)
  Up_JQ1 <- Up_JQ1[(Up_JQ1$Gene %in% Up) & (Up_JQ1$log2FoldChange < -0.5) & (Up_JQ1$padj < 0.05),]
  Up_JQ1 <- unlist(Up_JQ1$Gene)
  
  #Put into list
  gene_list <- list(Up, Down, Up_JQ1)
  names(gene_list) <- c(paste(agonist.name, 'Up', sep = "_"), 
                        paste(agonist.name, 'Down', sep = "_"), 
                        paste(agonist.name, 'Up_JQ1_Down', sep = "_"))
  return(gene_list)
}

# Get lists of all differentially upregulated/downregulated/JQ1 attenuated
PE_genes_up_down_jq1attenuated <- differential_genes(resIHW.PEonly, resIHW.pejq1.pedmso, 'PE')
ET1_genes_up_down_jq1attenuated <- differential_genes(resIHW.ET1only, resIHW.et1jq1.et1dmso, 'ET1')

# Venn diagrams of up and downregulated genes
# agonist.gene.set refers to the lists made with the differential_genes function
gene_venn_diagrams <- function(agonist.gene.set1, agonist.gene.set2, agonist.1, agonist.2, direction) {
                      # A - genes unique to agongist.gene.set1
  venn <- venneuler(c(A = length(setdiff(unlist(agonist.gene.set1[paste(agonist.1, direction, sep = "_")]), 
                                         unlist(agonist.gene.set2[paste(agonist.2, direction, sep = "_")]))),
                      # Genes that overlap between the two gene sets
                      'A&B' = length(intersect(unlist(agonist.gene.set1[paste(agonist.1, direction, sep = "_")]), 
                                               unlist(agonist.gene.set2[paste(agonist.2, direction, sep = "_")]))),
                      # B - genes that are unique to agonist.gene.set2
                      B = length(setdiff(unlist(agonist.gene.set2[paste(agonist.2, direction, sep = "_")]), 
                                         unlist(agonist.gene.set1[paste(agonist.1, direction, sep = "_")])))))
  plot(venn)
}

# Plot venn diagrams for each element of the gene.sets
venn_upregulated <- gene_venn_diagrams(PE_genes_up_down_jq1attenuated, ET1_genes_up_down_jq1attenuated, 'PE', 'ET1', 'Up')
venn_downregulated <- gene_venn_diagrams(PE_genes_up_down_jq1attenuated, ET1_genes_up_down_jq1attenuated, 'PE', 'ET1', 'Down')
venn_jq1_attenuated <- gene_venn_diagrams(PE_genes_up_down_jq1attenuated, ET1_genes_up_down_jq1attenuated, 'PE', 'ET1', 'Up_JQ1_Down')

# Heatmap of all genes abs(log2FC >1) and pvalue < 0.05
# Correct for batch effect with limma::removeBatchEffect
vsd <- vst(dds, blind=FALSE)
mat <- assay(vsd)
mat <- limma::removeBatchEffect(mat, vsd$batch)
colnames(mat) <- sample.table$sampleID
annoDF <- as.data.frame(colData(vsd)[,1,drop=FALSE])
rownames(annoDF) <- colnames(vsd)

# Extract only significant differentially expressed genes and Acetic.DMSO, PE.DMSO, and PE.JQ1 treatment groups
sig_genes_PE <- resIHW.PEonly[(abs(resIHW.PEonly$log2FoldChange) > 1) & (resIHW.PEonly$padj < 0.05),]
sig_genes_PE <- unlist(rownames(sig_genes_PE))
mat_PE <- mat[sig_genes_PE,]
mat_PE <- mat_PE[,c(1, 2, 9:12)]

set.seed(2)
# Find three K-means clusters
k_PE <- pheatmap(mat_PE, scale="row",kmeans_k = 3)

# Convert K-means clusters to a dataframe and rename column
clusterDF_PE <- as.data.frame(factor(k_PE$kmeans$cluster))
colnames(clusterDF_PE) <- "Cluster"

# Order the counts matrix by cluster number 
OrderByCluster_PE <- mat_PE[order(clusterDF_PE$Cluster),]

# Extract only significant differentially expressed genes and Acetic.DMSO, ET1.DMSO, and ET1.JQ1 treatment groups
sig_genes_ET1 <- resIHW.ET1only[(abs(resIHW.ET1only$log2FoldChange) > 1) & (resIHW.ET1only$padj < 0.05),]
sig_genes_ET1 <- unlist(rownames(sig_genes_ET1))
mat_ET1 <- mat[sig_genes_ET1,]
mat_ET1 <- mat_ET1[,c(1, 2, 5:8)]

set.seed(2)
# Find three K-means clusters
k_ET1 <- pheatmap(mat_ET1, scale="row", kmeans_k = 3)
# Convert K-means clusters to a dataframe and rename column
clusterDF_ET1 <- as.data.frame(factor(k_ET1$kmeans$cluster))
colnames(clusterDF_ET1) <- "Cluster"

#Renumber clusters to match PE clusters and sort data by them
clusterDF_ET1$Cluster <- ifelse(clusterDF_ET1$Cluster == 1, 4, clusterDF_ET1$Cluster)
clusterDF_ET1$Cluster <- ifelse(clusterDF_ET1$Cluster == 2, 1, clusterDF_ET1$Cluster)
clusterDF_ET1$Cluster <- ifelse(clusterDF_ET1$Cluster == 3, 2, clusterDF_ET1$Cluster)
clusterDF_ET1$Cluster <- ifelse(clusterDF_ET1$Cluster == 4, 3, clusterDF_ET1$Cluster)
clusterDF_ET1$Cluster <- factor(clusterDF_ET1$Cluster)
OrderByCluster_ET1 <- mat_ET1[order(clusterDF_ET1$Cluster),]

##Plot clustered heatmaps
PE_heatmap <- pheatmap(OrderByCluster_PE,
                       scale="row", annotation_row = clusterDF_PE,
                       show_rownames = FALSE,cluster_rows = FALSE, 
                       cluster_cols = TRUE, cellheight = 0.4)

ET1_heatmap <- pheatmap(OrderByCluster_ET1,
                        scale="row", annotation_row = clusterDF_ET1,
                        show_rownames = FALSE, cluster_rows = FALSE,
                        cluster_cols = TRUE, cellheight = 0.4)

# Combine the plots into one figure
p <- ggpubr::ggarrange(PE_heatmap[[4]], ET1_heatmap[[4]],
                       labels = c("A", "B"),
                       ncol = 2, nrow = 1)

# Save heatmaps
ggplot2::ggsave("pe_et1_heatmap_kmeans3_seed1.pdf", p, dpi = 600, width = 10, height = 6)





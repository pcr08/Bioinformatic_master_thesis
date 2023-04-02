#load packages
library(Seurat)
library(SingleCellSignalR)
require(Matrix)
require(stringr)
library(dplyr)

#load matrix
setwd('C:/Users/PAULA CAMARGO ROMERA/OneDrive - Uppsala universitet/Master/Uppsala University/Courses/Year 2/Thesis/processing')
data_newexp <- readRDS('data_newexp.rds')
dataset5 <- data_newexp@meta.data[data_newexp@meta.data$cell_type == "microglial cell",] # healthy microglial cell types
dataset6 <- data_newexp@meta.data[data_newexp@meta.data$cell_type == "macrophage",] # healthy macrophages cell types
dataset7 <- rbind(dataset5, dataset6) # put them together
names <- row.names(dataset7) #take the names of the genes that we want to keep
d1 <- subset(data_newexp, cells = names) #do a new subset with only interested cells
d1 <- SetIdent(d1, value = "cell_type")


## Quality Control (QC)
# plot histogram unique genes
nGenes <- d1$nFeature_RNA
#hist(nGenes, breaks = 50)
#abline(v = mean(nGenes), col = "black", lty = "longdash", lwd = 2)
# 2nd standard deviation
#abline(v = mean(nGenes) - 2*sd(nGenes), col = "red", lty = "longdash", lwd = 2)
#abline(v = mean(nGenes) + 2*sd(nGenes), col = "red", lty = "longdash", lwd = 2)
# plot histogram unique molecules
nMol <- d1$nCount_RNA
#hist(nMol, breaks = 50)
#abline(v = mean(nMol), col = "black", lty = "longdash", lwd = 2)
# 2nd standard deviation
#abline(v = mean(nMol) - 2*sd(nMol), col = "red", lty = "longdash", lwd = 2)
#abline(v = mean(nMol) + 2*sd(nMol), col = "red", lty = "longdash", lwd = 2)
#VlnPlot
#VlnPlot(d1, features = c("nFeature_RNA", "nCount_RNA"), group.by = "cell_type")
#FeatureScatter(object = d1, feature2 = 'nFeature_RNA', feature1 = 'nCount_RNA') #0.41
#Filtering
# Select cells (cell ids) with at least 5000 unique molecules
set1 <- d1$nCount_RNA >= 5000
set4 <- d1$nCount_RNA <= 30500
# Select cells (cell ids) with at least 200 unique genes
set2 <- d1$nFeature_RNA >= 200
# Select cells (cell ids) with as maximum 7000 unique genes (upper limit histogram). Upper limit is to make sure you don’t have any outliers that are dying cells or specific to just one or few samples. The cells or spots that are not uniformly distributed across all samples.
set3 <- d1$nFeature_RNA <= 7000
# Select cells where criteria above are met (set1 AND set2 AND set3 has to be TRUE)
keep.cells <- colnames(d1)[set1 & set2 & set3 & set4]
head(keep.cells)
# Now use these cell names to subset the "Seurat" object
d1.subset <- subset(d1, cells = keep.cells)
#Filter low expressed genes
# First get the expression matrix
sparse.exprMat <- GetAssayData(d1, slot = "counts")
# Calculate row sums
gene.counts <- Matrix::rowSums(sparse.exprMat)
# Check the results
head(gene.counts)
#expression gene higher than 1000
set5 <- gene.counts > 1000
head(set5)
keep.genes <- rownames(d1)[set5]
#subset
d1.subset <- subset(d1.subset, features = keep.genes) #123929 cells and 20385 genes (before 243377 cells in total and 28045 genes)
#FeatureScatter(object = d1.subset,feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA') #0.87


markers <- FindAllMarkers(d1.subset, log2FC.threshold = 0.2, test.use = "wilcox",
                          min.pct = 0.1, min.diff.pct = 0.2, only.pos = TRUE, max.cells.per.ident = 50,
                          assay = "RNA")

markers %>%
  group_by(cluster) %>%
  top_n(-5, p_val) -> top5_cell_selection

VlnPlot(d1.subset, features = as.character(unique(top5_cell_selection$gene)),
        ncol = 5, group.by = "cell_type", assay = "RNA", pt.size = 0.1)

View(markers)
write.csv(markers, file = "DEG.csv")


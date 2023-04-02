#load packages
library(Seurat)
library(SingleCellSignalR)
require(Matrix)
require(stringr)

#load matrix
setwd('C:/Users/PAULA CAMARGO ROMERA/OneDrive - Uppsala universitet/Master/Uppsala University/Courses/Year 2/Thesis/processing')
data_newexp <- readRDS('data_newexp.rds')
#four stages (MES/NPC/OPC/AC) of glioblastoma
dataset1 <- data_newexp@meta.data[data_newexp@meta.data$annotation_level_3 == "AC-like",]
dataset2 <- data_newexp@meta.data[data_newexp@meta.data$annotation_level_3 == "MES-like",]
dataset3 <- data_newexp@meta.data[data_newexp@meta.data$annotation_level_3 == "OPC-like",]
dataset4 <- data_newexp@meta.data[data_newexp@meta.data$annotation_level_3 == "NPC-like",]
dataset5 <- data_newexp@meta.data[data_newexp@meta.data$cell_type == "microglial cell",][1:25000,] # healthy microglial cell types
dataset6 <- data_newexp@meta.data[data_newexp@meta.data$cell_type == "macrophage",][1:35000,] # healthy macrophages cell types
dataset7 <- rbind(dataset1, dataset2) # put them together
dataset7 <- rbind(dataset7, dataset3) # put them together
dataset7 <- rbind(dataset7, dataset4) # put them together
dataset7 <- rbind(dataset7, dataset5) # put them together
dataset7 <- rbind(dataset7, dataset6) # put them together
names <- row.names(dataset7) #take the names of the genes that we want to keep
d1 <- subset(data_newexp, cells = names) #do a new subset with only interested cells
nCount = colSums(x = d1, slot = "counts")  # nCount_RNA
nFeature = colSums(x = GetAssayData(object = d1, slot = "counts") > 0)  # nFeatureRNA


## Quality Control (QC)
# plot histogram unique genes
nGenes <- d1$nFeature_RNA
hist(nGenes, breaks = 50)
abline(v = mean(nGenes), col = "black", lty = "longdash", lwd = 2)
# 2nd standard deviation
abline(v = mean(nGenes) - 2*sd(nGenes), col = "red", lty = "longdash", lwd = 2)
abline(v = mean(nGenes) + 2*sd(nGenes), col = "red", lty = "longdash", lwd = 2)
# plot histogram unique molecules
nMol <- d1$nCount_RNA
hist(nMol, breaks = 50)
abline(v = mean(nMol), col = "black", lty = "longdash", lwd = 2)
# 2nd standard deviation
abline(v = mean(nMol) - 2*sd(nMol), col = "red", lty = "longdash", lwd = 2)
abline(v = mean(nMol) + 2*sd(nMol), col = "red", lty = "longdash", lwd = 2)
#VlnPlot
VlnPlot(d1, features = c("nFeature_RNA", "nCount_RNA"), group.by = "cell_type")
FeatureScatter(object = d1, feature2 = 'nFeature_RNA', feature1 = 'nCount_RNA') #0.41
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
d1.subset <- subset(d1.subset, features = keep.genes)
FeatureScatter(object = d1.subset,feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA')

#normalization
library(glmGamPoi)
d1.subset <- SCTransform(object = d1.subset, vars.to.regress = "cell_type", method="glmGamPoi", conserve.memory = TRUE) #includes scaling and transformation in one step

#variable feature identification
d1.subset <- FindVariableFeatures(d1.subset, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
top10 <- head(VariableFeatures(d1.subset), 10)
plot1 <- VariableFeaturePlot(d1.subset)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE) #+ theme(legend.position="top")
plot2

#dimensionality reduction
d1.subset <- ScaleData(d1.subset, vars.to.regress = "cell_type")
d1.subset <- RunPCA(d1.subset, features = VariableFeatures(object = d1.subset), verbose = FALSE)
ElbowPlot(d1.subset) #to see how many PCA explains the variation -- around 15
VizDimLoadings(d1.subset, dims = 1:2, reduction = "pca")
DimPlot(d1.subset, reduction = "pca", group.by = c("cell_type","annotation_level_3"))
DimHeatmap(d1.subset, dims = 1:14, cells = 500, balanced = TRUE)

#non-linear dimensional reduction
#d1.subset <- RunUMAP(d1.subset, dims = 1:30, reduction = "pca")
#DimPlot(d1.subset, reduction = "umap", group.by = "replicate", label = TRUE)
#d1.subset <- RunTSNE(d1.subset, dims = 1:20, verbose = FALSE)
#DimPlot(pbmc, reduction = "tsne", label = TRUE)


#translate Ensembl symbols to HUGO
library("EnsDb.Hsapiens.v86")
mmatrix <- GetAssayData(object = d1.subset, assay= "SCT", slot = "data") #get the normalized UMI counts for both GBM and microglial cells
row.names(mmatrix) <- mapIds(EnsDb.Hsapiens.v86,keys = row.names(mmatrix), keytype = "GENEID",column = "SYMBOL", multiVals="first")

which(is.na(rownames(mmatrix)))
mmatrix <- mmatrix[-c(19155,19180,19505,19811,19811,19812,19816,19822,19831,19833,19839,19841,19842,19846,19850,19853,19855,19857,19861,19867),,drop=FALSE]


#singlecellsignalR
library(SCOPfunctions)
data <- utils_big_as.matrix(mmatrix)
clust <- as.numeric(d1.subset$annotation_level_3)
clust[clust == 1] = 3 #glioblastoma AC-like
clust[clust == 20] = 1 #microglial cell (TAM-MG)
clust[clust == 19] = 2 #macrophages (TAM-BDM)
clust[clust == 7] = 4 #glioblastoma MES-like
clust[clust == 12] = 5 #glioblastoma NPC-like
clust[clust == 15] = 6 #glioblastoma OPC-like
clust.ana <- cluster_analysis(data = data, c.names = c("microglial cells", "macrophages", "glioblastoma AC-like", "glioblastoma MES-like", "glioblastoma NPC-like", "glioblastoma OPC-like"), genes = rownames(data), cluster = clust)
signal <- cell_signaling(data = data, c.names = c("microglial cells", "macrophages", "glioblastoma AC-like", "glioblastoma MES-like", "glioblastoma NPC-like", "glioblastoma OPC-like"), genes = rownames(data), cluster = clust)
#install.packages("igraph")
library(igraph) #for modified function
environment(inter_network) <- asNamespace('SingleCellSignalR') #modified function
assignInNamespace("inter_network", inter_network, ns = "SingleCellSignalR")
inter.net <- inter_network(data = data, c.names = c("microglial cells", "macrophages", "glioblastoma AC-like", "glioblastoma MES-like", "glioblastoma NPC-like", "glioblastoma OPC-like"), signal = signal, genes = rownames(data), cluster = clust)
SingleCellSignalR::visualize_interactions(signal, write.out = TRUE)

#load packages
library(Seurat)
library(SingleCellSignalR)
require(Matrix)
require(stringr)

#load matrix
setwd('/proj/snic2022-22-1049/Paula/GBmap')
data_newexp <- readRDS('data_newexp.rds')
setwd('/proj/snic2022-22-1049/Paula/GBmap/othercelltypes')
#glioblastoma vs other types
dataset1 <- data_newexp@meta.data[data_newexp@meta.data$cell_type == "malignant cell",]
dataset2 <- data_newexp@meta.data[data_newexp@meta.data$cell_type == "mast cell",]
dataset3 <- data_newexp@meta.data[data_newexp@meta.data$cell_type == "endothelial cell",]
dataset4 <- data_newexp@meta.data[data_newexp@meta.data$cell_type == "astrocyte",]
dataset5 <- data_newexp@meta.data[data_newexp@meta.data$cell_type == "oligodendrocyte",]
dataset6 <- data_newexp@meta.data[data_newexp@meta.data$cell_type == "dendritic cell",]
dataset7 <- data_newexp@meta.data[data_newexp@meta.data$cell_type == "neuron",]
dataset8 <- data_newexp@meta.data[data_newexp@meta.data$cell_type == "monocyte",]
dataset9 <- data_newexp@meta.data[data_newexp@meta.data$cell_type == "natural killer cell",]
dataset10 <- data_newexp@meta.data[data_newexp@meta.data$cell_type == "radial glial cell",]
dataset11 <- data_newexp@meta.data[data_newexp@meta.data$cell_type == "mural cell",]

dataset12 <- rbind(dataset2, dataset3) # put them together
dataset12 <- rbind(dataset12, dataset4) # put them together
dataset12 <- rbind(dataset12, dataset5) # put them together
dataset12 <- rbind(dataset12, dataset6) # put them together
dataset12 <- rbind(dataset12, dataset7) # put them together
dataset12 <- rbind(dataset12, dataset8) # put them together
dataset12 <- rbind(dataset12, dataset9) # put them together
dataset12 <- rbind(dataset12, dataset10) # put them together
dataset12 <- rbind(dataset12, dataset11) # put them together
names <- row.names(dataset1) #take the names of the genes that we want to keep, malignant cells
d1 <- subset(data_newexp, cells = names) #do a new subset with only interested cells
nCount = colSums(x = d1, slot = "counts")  # nCount_RNA
nFeature = colSums(x = GetAssayData(object = d1, slot = "counts") > 0)  # nFeatureRNA
names2 <- row.names(dataset12) #other cell types
d2 <- subset(data_newexp, cells = names2)

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


se_combined.subset <- merge(d1.subset, y = d2)


#normalization
library(glmGamPoi)
se_combined.subset <- SCTransform(object = se_combined.subset, vars.to.regress = "cell_type", method="glmGamPoi", conserve.memory = TRUE) #includes scaling and transformation in one step

#variable feature identification
se_combined.subset <- FindVariableFeatures(se_combined.subset, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
#dimensionality reduction
se_combined.subset <- ScaleData(se_combined.subset, vars.to.regress = "cell_type")
se_combined.subset <- RunPCA(se_combined.subset, features = VariableFeatures(object = se_combined.subset), verbose = FALSE)

mmatrix <- GetAssayData(object = se_combined.subset, assay= "SCT", slot = "data") #get the normalized UMI counts for both GBM and microglial cells

library("EnsDb.Hsapiens.v86")
row.names(mmatrix) <- mapIds(EnsDb.Hsapiens.v86,keys = row.names(mmatrix), keytype = "GENEID",column = "SYMBOL", multiVals="first")

which(is.na(rownames(mmatrix)))
mmatrix <- mmatrix[-c(18618,18643,18965,19261,19262,19266,19270,19283,19286,19287,19309,19311,19313,19314,19315,19318,19324,19326,19340,19343,19344,19345,19346,19347,19348,19363,19350,19351,19354,19355,19357,19361,19363,19365,19369,19376,25213,25228,25249,25251,25254,25265,25272,25273,25277,25282,25285,25286),,drop=FALSE]
which(is.na(rownames(mmatrix)))

#singlecellsignalR
library(SCOPfunctions)
data <- utils_big_as.matrix(mmatrix)
x1<-rep(c(1),each=61408) #malignant cell 61408
x2<-rep(c(2),each=373) #mast cell 373
x3<-rep(c(3),each=673) #endothelial cell 673
x4<-rep(c(4),each=173) #astrocyte 173
x5<-rep(c(5),each=12481) #oligodendrocyte 12481
x6<-rep(c(6),each=3961) #dendritic cell 3961
x7<-rep(c(7),each=22) #neuron 22
x8<-rep(c(8),each=14215) #monocyte 14215
x9<-rep(c(9),each=2489) #natural killer cell 2489
x10<-rep(c(10),each=2807) #radial glial cell 2807
x11<-rep(c(11),each=1418) #mural cell 1418
clust <- c(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11)
#clust.ana <- cluster_analysis(data = data, c.names = c("glioblastoma","mast cells", "endothelial cells","astrocytes", "oligodendrocytes", "dendritic cells", "neurons", "monocytes", "natural killer cells", "radial glial cell","mural cells"), genes = rownames(data), cluster = clust)
signal <- cell_signaling(data = data, c.names = c("glioblastoma","mast cells", "endothelial cells","astrocytes", "oligodendrocytes", "dendritic cells", "neurons", "monocytes", "natural killer cells", "radial glial cell","mural cells"), genes = rownames(data), cluster = clust, s.score=0)
#install.packages("igraph")
library(igraph) #for modified function
environment(inter_network) <- asNamespace('SingleCellSignalR') #modified function
assignInNamespace("inter_network", inter_network, ns = "SingleCellSignalR")
inter.net <- inter_network(data = data, c.names = c("glioblastoma","mast cells", "endothelial cells","astrocytes", "oligodendrocytes", "dendritic cells", "neurons", "monocytes", "natural killer cells", "radial glial cell","mural cells"), signal = signal, genes = rownames(data), cluster = clust)
#SingleCellSignalR::visualize_interactions(signal, write.out = TRUE)



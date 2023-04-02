library(boot)

function_data <- function(data, i){
  d2 <- data[i,]
  list_d <- tail(n=50, mmatrix@Dimnames[[2]])
  j=26166
  for (p in 1:50){
    if (p < 50){
      if (p == 1){
        l <- 0.9*d2[,j] + 0.1*d2[,j+1]
      }else{
        z <- 0.9*d2[,j] + 0.1*d2[,j+1]
        l <- rbind(l,z)
      }
    }else{
      z <- 0.9*d2[,j] + 0.1*d2[26166] #last value, compare it with first
      l <- rbind(l,z)
    }
    j <- j + 1
  }
  return(l)
}

setwd('C:/Users/PAULA CAMARGO ROMERA/OneDrive - Uppsala universitet/Master/Uppsala University/Courses/Year 2/Thesis/processing/our_data')
dataset1 <- read.table(file = "table_paula.txt", header = TRUE, row.names = NULL)
#set.seed(1)
#bootstrap_1 <- boot(dataset1,function_1,R=500) #let's do 500 repeats so it is the same amount of microglia than glioblastoma (26165 vs 25000, 50x500 = 25000)
#bootstrap <- boot(df,function_2,R=500)

#load packages
library(Seurat)
library(SingleCellSignalR)
require(Matrix)
require(stringr)

#load matrix
setwd('C:/Users/PAULA CAMARGO ROMERA/OneDrive - Uppsala universitet/Master/Uppsala University/Courses/Year 2/Thesis/processing/our_data')
dataset1 <- read.table(file = "table_paula.txt", header = TRUE, row.names = NULL)
rownames(dataset1) <- make.names(dataset1[,1], unique =  T)
dataset1[,1] <- NULL
first_column <- colnames(dataset1)
second_column <- rep(c("neoplastic"),each=length(colnames(dataset1)))
#log2(x/sum(x)*10000), ie if xi is column i of x, then log2(1 +   xi/sum(xi)*10000)
metadata1 <- data.frame(first_column, second_column)
colnames(metadata1) <- c('sample','cell_type')
rownames(metadata1) <- make.names(metadata1[,1], unique =  T)
se_neo <- CreateSeuratObject(counts = dataset1, meta.data = metadata1) #create Seurat object for lab data (tumour)
se_neo <- NormalizeData(se_neo, normalization.method = "LogNormalize" )


setwd('C:/Users/PAULA CAMARGO ROMERA/OneDrive - Uppsala universitet/Master/Uppsala University/Courses/Year 2/Thesis/processing')
data_newexp <- readRDS('data_newexp.rds')
dataset2 <- data_newexp@meta.data[data_newexp@meta.data$cell_type == "microglial cell",] # healthy microglial cell types
se_m <- subset(data_newexp, cells = row.names(dataset2))
#nCount = colSums(x = se_m, slot = "counts")  # nCount_RNA
#nFeature = colSums(x = GetAssayData(object = se_m, slot = "counts") > 0)  # nFeatureRNA

#QC for each dataset
#our data -- neoplastic cells
# plot histogram unique genes
nGenes <- se_neo$nFeature_RNA
#hist(nGenes, breaks = 50)
#abline(v = mean(nGenes), col = "black", lty = "longdash", lwd = 2)
# 2nd standard deviation
#abline(v = mean(nGenes) - 2*sd(nGenes), col = "red", lty = "longdash", lwd = 2)
#abline(v = mean(nGenes) + 2*sd(nGenes), col = "red", lty = "longdash", lwd = 2)
# plot histogram unique molecules
nMol <- se_neo$nCount_RNA
#hist(nMol, breaks = 50)
#abline(v = mean(nMol), col = "black", lty = "longdash", lwd = 2)
# 2nd standard deviation
#abline(v = mean(nMol) - 2*sd(nMol), col = "red", lty = "longdash", lwd = 2)
#abline(v = mean(nMol) + 2*sd(nMol), col = "red", lty = "longdash", lwd = 2)
#VlnPlot
#VlnPlot(se_neo, features = c("nFeature_RNA", "nCount_RNA"), group.by = "cell_type")
#FeatureScatter(object = se_neo, feature2 = 'nFeature_RNA', feature1 = 'nCount_RNA') #0.56, I don't filter because I don't want to lose patients from the lab
sparse.exprMat <- GetAssayData(se_neo, slot = "counts")
# Calculate row sums
gene.counts_neo <- Matrix::rowSums(sparse.exprMat)
# Check the results
head(gene.counts_neo)
#expression gene higher than 500
set5 <- gene.counts_neo > 1000
head(set5)
keep.genes_neo <- rownames(se_neo)[set5]
#subset
se_neo.subset <- subset(se_neo, features = keep.genes_neo)
#FeatureScatter(object = se_neo.subset,feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA') 





#GBmap data -- microglial cells
# plot histogram unique genes
nGenes <- se_m$nFeature_RNA
#hist(nGenes, breaks = 50)
#abline(v = mean(nGenes), col = "black", lty = "longdash", lwd = 2)
# 2nd standard deviation
#abline(v = mean(nGenes) - 2*sd(nGenes), col = "red", lty = "longdash", lwd = 2)
#abline(v = mean(nGenes) + 2*sd(nGenes), col = "red", lty = "longdash", lwd = 2)
# plot histogram unique molecules
nMol <- se_m$nCount_RNA
#hist(nMol, breaks = 50)
#abline(v = mean(nMol), col = "black", lty = "longdash", lwd = 2)
# 2nd standard deviation
#abline(v = mean(nMol) - 2*sd(nMol), col = "red", lty = "longdash", lwd = 2)
#abline(v = mean(nMol) + 2*sd(nMol), col = "red", lty = "longdash", lwd = 2)
#VlnPlot
#VlnPlot(se_m, features = c("nFeature_RNA", "nCount_RNA"), group.by = "cell_type")
#FeatureScatter(object = se_m, feature2 = 'nFeature_RNA', feature1 = 'nCount_RNA') #0.23
#Filtering
#convert it to bulk grouping by samples (only taking the microglial) and then perform the analysis with the tumour
# Select cells (cell ids) with at least 5000 unique molecules of microglial cell
set1 <- se_m$nCount_RNA >= 5000
set4 <- se_m$nCount_RNA <= 30500
# Select cells (cell ids) with at least 200 unique genes of microglial cell
set2 <- se_m$nFeature_RNA >= 200
# Select cells (cell ids) with as maximum 7000 unique genes (upper limit histogram) of microglial cell. Upper limit is to make sure you don’t have any outliers that are dying cells or specific to just one or few samples. The cells or spots that are not uniformly distributed across all samples.
set3 <- se_m$nFeature_RNA <= 7000
# Select cells where criteria above are met (set1 AND set2 AND set3 AND set4 for microgial cell, OR set6 AND set7 for tumour cell)
keep.cells <- colnames(se_m)[(set1 & set2 & set3 & set4)]
head(keep.cells)
# Now use these cell names to subset the "Seurat" object
se_m.subset <- subset(se_m, cells = keep.cells)
#Filter low expressed genes
# First get the expression matrix
sparse.exprMat <- GetAssayData(se_m, slot = "data")
# Calculate row sums
gene.counts <- Matrix::rowSums(sparse.exprMat)
# Check the results
head(gene.counts)
#expression gene higher than 1000
set5 <- gene.counts > 1000
head(set5)
keep.genes <- rownames(se_m)[set5]
#subset
se_m.subset <- subset(se_m.subset, features = keep.genes)
#FeatureScatter(object = se_m.subset,feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA')


#combine datasets
se_combined.subset <- merge(se_m.subset, y = se_neo.subset, add.cell.ids = c("microglial cell", "neoplastic"))

#normalization
library(glmGamPoi)
se_combined.subset <- SCTransform(object = se_combined.subset, vars.to.regress = "cell_type", method="glmGamPoi", conserve.memory = TRUE) #includes scaling and transformation in one step

#variable feature identification
se_combined.subset <- FindVariableFeatures(se_combined.subset, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

#dimensionality reduction
se_combined.subset <- ScaleData(se_combined.subset, vars.to.regress = "cell_type")
se_combined.subset <- RunPCA(se_combined.subset, features = VariableFeatures(object = se_combined.subset), verbose = FALSE)
#d1.subset <- RunPCA(d1.subset, verbose = FALSE)
#ElbowPlot(se_combined.subset) #to see how many PCA explains the variation -- around 15
#VizDimLoadings(se_combined.subset, dims = 1:2, reduction = "pca")
#DimPlot(se_combined.subset, reduction = "pca", group.by = "cell_type")
#DimHeatmap(se_combined.subset, dims = 1:14, cells = 500, balanced = TRUE)


#translate Ensembl symbols to HUGO
library("EnsDb.Hsapiens.v86")
#library("org.Hs.eg.db")
mmatrix <- GetAssayData(object = se_combined.subset, assay= "SCT", slot = "data") #get the normalized UMI counts for both GBM and microglial cells
m1 <- row.names(mmatrix[1:8997,]) #microglial, ENSEMBL
m2 <- rownames(mmatrix[8998:24489,]) #lab, HUGO
m1 <- mapIds(EnsDb.Hsapiens.v86,keys = m1, keytype = "GENEID",column = "SYMBOL", multiVals="first")
which(is.na(rownames(m1)))
m3 <- c(m1,m2)
row.names(mmatrix) <- m3
which(is.na(rownames(mmatrix)))

#bootstrap1 <- boot(data = data,function_data(data = data, j = 24439),R=1)

#singlecellsignalR
library(SCOPfunctions)
data <- utils_big_as.matrix(mmatrix)

bootstrap1 <- boot(data,statistic = function_data,R=500) #repeat 500 times, so we get 25000 columns for glioblastoma patients
for (s in 1:500){
  j <- 1
  for (i in 1:50){ #split the final row to 50 patients
    if (i != 50){
      crow <- bootstrap1$t[s,j:(j+24488)]
      data <- cbind(data, crow)
      j <- j + 24489
    }else{
      crow <- bootstrap1$t[s,j:1224450]
      data <- cbind(data, crow)
    }
  }
}

names1 <- colnames(data)
names2 <- names1[1:26215]
rnames <- c(1:25000)
names3 <- c(names2, rnames)
colnames(data) <- names3

clust <- se_combined.subset$cell_type
clust[clust == "neoplastic"] = 1 #glioblastoma 50
clust[clust == "microglial cell"] = 2 #microglial cell 26165
x1<-rep(c(1),each=25000) #add new glioblastoma values
clust <- append(clust,x1)
clust.ana <- cluster_analysis(data = data, c.names = c("glioblastoma","microglial cells"), genes = rownames(data), cluster = clust)
signal <- cell_signaling(data = data, c.names = c("glioblastoma","microglial cells"), genes = rownames(data), cluster = clust)
#install.packages("igraph")
library(igraph) #for modified function
environment(inter_network) <- asNamespace('SingleCellSignalR') #modified function
assignInNamespace("inter_network", inter_network, ns = "SingleCellSignalR")
inter.net <- inter_network(data = data, c.names = c("glioblastoma","microglial cells"), signal = signal, genes = rownames(data), cluster = clust)
SingleCellSignalR::visualize_interactions(signal, write.out = TRUE)


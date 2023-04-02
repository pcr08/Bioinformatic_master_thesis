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
dataset2 <- data_newexp@meta.data[data_newexp@meta.data$cell_type == "dendritic cell",] # dendritic cell types
se_m <- subset(data_newexp, cells = row.names(dataset2))

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
FeatureScatter(object = se_neo, feature2 = 'nFeature_RNA', feature1 = 'nCount_RNA') #0.56, I don't filter because I don't want to lose patients from the lab
sparse.exprMat <- GetAssayData(se_neo, slot = "counts")
# Calculate row sums
gene.counts_neo <- Matrix::rowSums(sparse.exprMat)
# Check the results
head(gene.counts_neo)
#expression gene higher than 1000
set5 <- gene.counts_neo > 1000
head(set5)
keep.genes_neo <- rownames(se_neo)[set5]
#subset
se_neo.subset <- subset(se_neo, features = keep.genes_neo)
FeatureScatter(object = se_neo.subset,feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA') 


#GBmap data -- dendritic cells
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
#FeatureScatter(object = se_m, feature2 = 'nFeature_RNA', feature1 = 'nCount_RNA')
#Filtering
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
se_m.subset <- subset(se_m, features = keep.genes)
#FeatureScatter(object = se_m.subset,feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA')


#combine datasets
se_combined.subset <- merge(se_m.subset, y = se_neo.subset, add.cell.ids = c("dendritic_cells", "neoplastic"))


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
mmatrix <- GetAssayData(object = se_combined.subset, assay= "SCT", slot = "data") #get the normalized UMI counts for both GBM and dendritic cells


library(SCOPfunctions)
data <- utils_big_as.matrix(mmatrix)

nc <- ncol(data[,1:3961]) + 3962
new_mat <- matrix(0, ncol=nc,nrow=nrow(data[,1:3961]))
new_mat[1:nrow(data), 1:ncol(data[,1:3961])] <- data[,1:3961] #new matrix with the values from old one (with dendritic cells) and the rest 0 to add noise

number_fake_cells = 3961
gamma = 0.9
r1 <- 3963
c1 <- 7923
list_d <- tail(n=50, mmatrix@Dimnames[[2]])
p = 0
setwd('C:/Users/PAULA CAMARGO ROMERA/OneDrive - Uppsala universitet/Master/Uppsala University/Courses/Year 2/Thesis/processing/our_data')
for (patient in 3962:4011){
  Sys.time()
  p = p + 1
  fake_cells <- matrix(0, nrow=length(data[,1]), ncol=number_fake_cells) #new matrix with noise to be added
  for (i in 1:number_fake_cells){
    rpatient <- sample(1:ncol(data), 1, replace=TRUE)
    fake_cells[,i] <- gamma*data[,patient]+(1-gamma)*data[,rpatient] #vector with new values
  }
  new_mat[,3962] <- data[,patient] #add patient
  new_mat[,r1:c1] <- fake_cells #add noise to the original matrix
  
  names1 <- colnames(data[,1:3961]) #dendritic cells
  names2 <- colnames(data[,(patient-1):patient]) #patient from glioblastoma
  names2 <- names2[2]
  rnames <- c(1:3961) #glioblastoma
  names3 <- c(names1, names2)
  names4 <- c(names3, rnames)
  colnames(new_mat) <- names4
  rownames(new_mat) <- rownames(data[,1:2])
  Sys.time()
  folder <- paste('/proj/snic2022-22-1049/Paula/new_boot_CPDB/dendritic_cells/', list_d[[p]], sep="")
  dir.create(folder)
  setwd(folder)
  file <- paste('matrix_cellphonedb_normalisation_othercelltypes_ind_', list_d[[p]], sep="")
  file1 <- paste(file, '.txt', sep="")
  write.table(new_mat, file1, append = FALSE, sep = "\t", dec = ".",
              row.names = TRUE, col.names = TRUE)
}

# #metadata normalisation
# setwd('C:/Users/PAULA CAMARGO ROMERA/OneDrive - Uppsala universitet/Master/Uppsala University/Courses/Year 2/Thesis/processing/cellphonedb_results/gbmap_celltypes/dendritic_cells')
# p <- 0
# for (j in 3962:4011){
#   p <- p + 1
#   x1<-rep(c("dendritic_cells"),each=3961) #dendritic cells
#   x2<-rep(c("neoplastic"),each=3962) #glioblastoma 1 patients lab data + noise
# 
#   l1 <- str_split(se_combined.subset@assays$SCT@data@Dimnames[[2]], '"')
#   dd1 <- c(l1[1:3961])
#   dd2 <- c(l1[j])
#   v1 <- c(1:3961) #random names for the noise
#   dd2 <- c(dd2, v1)
# 
#   df1 <- data.frame(unlist(dd1),x1)
#   df2 <- data.frame(unlist(dd2),x2)
#   names(df1) <- c('cell', 'cell_type')
#   names(df2) <- c('cell', 'cell_type')
#   df3 <- rbind(df1, df2) # put them together
#   file <- paste('meta_cellphone_normalisation_ourdata_bootstrap_ind', list_d[[p]], sep="")
#   file1 <- paste(file, ".txt", sep="")
#   write.table(df3,file1, row.names = FALSE, sep = "\t")
# }

#post commands -- python + linux commands

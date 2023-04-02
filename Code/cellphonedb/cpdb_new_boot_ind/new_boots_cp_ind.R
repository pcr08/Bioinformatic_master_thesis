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
#expression gene higher than 500
set5 <- gene.counts_neo > 1000
head(set5)
keep.genes_neo <- rownames(se_neo)[set5]
#subset
se_neo.subset <- subset(se_neo, features = keep.genes_neo)
FeatureScatter(object = se_neo.subset,feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA') 


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
se_combined.subset <- merge(se_m.subset, y = se_neo.subset, add.cell.ids = c("microglial_cell", "neoplastic"))


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


#mmatrix <- GetAssayData(object = d1, slot = "data")[, names] #get the normalized UMI counts for both GBM and microglial cells
#translate Ensembl symbols to HUGO
library("EnsDb.Hsapiens.v86")
#library("org.Hs.eg.db")
mmatrix <- GetAssayData(object = se_combined.subset, assay= "SCT", slot = "data") #get the normalized UMI counts for both GBM and microglial cells


library(SCOPfunctions)
data <- utils_big_as.matrix(mmatrix)
#data1 <- data[,1:15000] #take only microglial cells, not glioblastoma patients

nc <- ncol(data[,1:15000]) + 15001
new_mat <- matrix(0, ncol=nc,nrow=nrow(data[,1:15000]))
new_mat[1:nrow(data), 1:ncol(data[,1:15000])] <- data[,1:15000] #new matrix with the values from old one (will microglia) and the rest 0 to add noise

number_fake_cells = 15000
gamma = 0.9
r1 <- 15002
c1 <- 30001
list_d <- tail(n=50, mmatrix@Dimnames[[2]])
p = 34
setwd('C:/Users/PAULA CAMARGO ROMERA/OneDrive - Uppsala universitet/Master/Uppsala University/Courses/Year 2/Thesis/processing/our_data')
#patients has to be from 26166:26215
for (patient in 26200:26215){
  Sys.time()
  p = p + 1
  fake_cells <- matrix(0, nrow=length(data[,1]), ncol=number_fake_cells) #new matrix with noise to be added
  for (i in 1:number_fake_cells){
    #rpatient <- runif(1, 1, ncol(m1))
    rpatient <- sample(1:ncol(data), 1, replace=TRUE)
    fake_cells[,i] <- gamma*data[,patient]+(1-gamma)*data[,rpatient] #vector with new values
  }
  new_mat[,15001] <- data[,patient] #add patient
  new_mat[,r1:c1] <- fake_cells #add noise to the original matrix
  
  names1 <- colnames(data[,1:15000]) #microglial
  names2 <- colnames(data[,(patient-1):patient]) #patient from glioblastoma
  names2 <- names2[2]
  rnames <- c(1:15000) #glioblastoma
  names3 <- c(names1, names2)
  names4 <- c(names3, rnames)
  colnames(new_mat) <- names4
  rownames(new_mat) <- rownames(data[,1:2])
  Sys.time()
  file <- paste('matrix_cellphonedb_normalisation_ourdata_bootstrap_ind', list_d[[p]], sep="")
  file1 <- paste(file, '.txt', sep="")
  write.table(new_mat, file1, append = FALSE, sep = "\t", dec = ".",
              row.names = TRUE, col.names = TRUE)
}

# #metadata normalisation
# setwd('c:/users/paula camargo romera/onedrive - uppsala universitet/master/uppsala university/courses/year 2/thesis/processing/our_data/cellphonedb_results/decrease')
# p <- 0
# for (j in 26166:26215){
#   p <- p + 1
#   x1<-rep(c("microglial cell"),each=15000) #microglial cell (tam-mg) 26165
#   x2<-rep(c("neoplastic"),each=15001) #glioblastoma 1 patients lab data + 25000 noise
# 
#   l1 <- str_split(se_combined.subset@assays$SCT@data@Dimnames[[2]], '"')
#   dd1 <- c(l1[1:15000])
#   dd2 <- c(l1[j])
#   v1 <- c(1:15000) #random names for the noise
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

#post commands
#python commands
#with open('C:/Users/PAULA CAMARGO ROMERA/OneDrive - Uppsala universitet/Master/Uppsala University/Courses/Year 2/Thesis/processing/our_data/cellphonedb_results/bootstrap/meta_cellphone_normalisation_ourdata_bootstrap_indneoplastic_U3299.txt', 'r', encoding='utf-8') as input_file:
#  lines = input_file.readlines()
#print(lines)
#with open('C:/Users/PAULA CAMARGO ROMERA/OneDrive - Uppsala universitet/Master/Uppsala University/Courses/Year 2/Thesis/processing/our_data/cellphonedb_results/bootstrap/meta_cellphone2_normalisation_ourdata_bootstrap_indneoplastic_U3299.txt', 'w', encoding='utf-8') as output_file:
#  for line in lines:
#  output_file.write(line.replace('microglial\tcell', 'microglial_cell'))

#linux commands wsl
#for file in 'C:/Users/PAULA CAMARGO ROMERA/OneDrive - Uppsala universitet/Master/Uppsala University/Courses/Year 2/Thesis/processing/our_data/cellphonedb_results/bootstrap/decrease/meta/*'
#do
#f1="${file}2.txt"
#sed -e 's/"//g' file > f1
#done

#sed -e 's/"//g' meta_cellp_normalisation_ourdata_bootstrap.txt > meta_cellp_normalisation_ourdata_bootstrap2.txt


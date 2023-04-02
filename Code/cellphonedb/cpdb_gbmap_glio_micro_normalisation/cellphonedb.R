library(SingleCellSignalR)
require(Seurat)
require(Matrix)
require(stringr)
#cellphonedb with SCdata_steps normalization
#load matrix
setwd('C:/Users/PAULA CAMARGO ROMERA/OneDrive - Uppsala universitet/Master/Uppsala University/Courses/Year 2/Thesis/processing')
data_newexp <- readRDS('data_newexp.rds')
dataset1 <- data_newexp@meta.data[data_newexp@meta.data$annotation_level_1 == "Neoplastic",] # to get the ones that are GBM
dataset2 <- data_newexp@meta.data[data_newexp@meta.data$cell_type == "microglial cell",] # healthy microglial cell types
dataset3 <- rbind(dataset1, dataset2) # put them together 
names <- row.names(dataset3) #take the names of the genes that we want to keep
d1 <- subset(data_newexp, cells = names) #do a new subset with only interested cells
nCount = colSums(x = d1, slot = "counts")  # nCount_RNA
nFeature = colSums(x = GetAssayData(object = d1, slot = "counts") > 0)  # nFeatureRNA
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
#expression gene higher than 500
set5 <- gene.counts > 1000
head(set5)
keep.genes <- rownames(d1)[set5]
#subset
d1.subset <- subset(d1.subset, features = keep.genes)
#normalization
library(glmGamPoi)
d1.subset <- SCTransform(object = d1.subset, vars.to.regress = "cell_type", method="glmGamPoi", conserve.memory = TRUE) #includes scaling and transformation in one step
#variable feature identification
d1.subset <- FindVariableFeatures(d1.subset, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
#dimensionality reduction
d1.subset <- ScaleData(d1.subset, vars.to.regress = "cell_type")
d1.subset <- RunPCA(d1.subset, features = VariableFeatures(object = d1.subset), verbose = FALSE)
mmatrix <- GetAssayData(object = d1.subset, assay= "SCT", slot = "data") #get the normalized UMI counts for both GBM and microglial cells
library(SCOPfunctions)
data <- utils_big_as.matrix(mmatrix)
file <- "matrix_cellphonedb_normalisation.txt"
write.table(data, file, append = FALSE, sep = "\t", dec = ".",
            row.names = TRUE, col.names = TRUE)

#metadata normalisation
setwd('C:/Users/PAULA CAMARGO ROMERA/OneDrive - Uppsala universitet/Master/Uppsala University/Courses/Year 2/Thesis/processing')
x1<-rep(c("Neoplastic"),each=61408) #sum(d1.subset@meta.data$cell_type == "malignant cell")
x2<-rep(c("microglial cell"),each=26165) #sum(d1.subset@meta.data$cell_type == "microglial cell")
l1 <- str_split(d1.subset@assays$SCT@data@Dimnames[[2]], '"')
dd1 <- c(l1[1:61408])
dd2 <- c(l1[61409:87573])
df1 <- data.frame(unlist(dd1),x1)
df2 <- data.frame(unlist(dd2),x2)
names(df1) <- c('Cell', 'cell_type')
names(df2) <- c('Cell', 'cell_type')
df3 <- rbind(df1, df2) # put them together
write.table(df3,"meta_cellphone_normalisation.txt", row.names = FALSE, sep = "\t")



#post commands
#python commands
#with open('C:/Users/PAULA CAMARGO ROMERA/OneDrive - Uppsala universitet/Master/Uppsala University/Courses/Year 2/Thesis/processing/meta_cellphone_normalisation.txt.txt', 'r', encoding='utf-8') as input_file:
#  lines = input_file.readlines()
#print(lines)
#with open('C:/Users/PAULA CAMARGO ROMERA/OneDrive - Uppsala universitet/Master/Uppsala University/Courses/Year 2/Thesis/processing/meta_normalisation.txt', 'w', encoding='utf-8') as output_file:
#  for line in lines:
#  output_file.write(line.replace('microglial\tcell', 'microglial_cell'))
#linux commands wsl
#sed -e 's/"//g' meta_normalisation.txt > meta_normalisation2.txt

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
dataset5 <- data_newexp@meta.data[data_newexp@meta.data$cell_type == "microglial cell",] # healthy microglial cell types
dataset6 <- data_newexp@meta.data[data_newexp@meta.data$cell_type == "macrophage",] # healthy macrophages cell types
dataset7 <- rbind(dataset1, dataset2) # put them together
dataset7 <- rbind(dataset7, dataset3) # put them together
dataset7 <- rbind(dataset7, dataset4) # put them together
dataset7 <- rbind(dataset7, dataset5) # put them together
dataset7 <- rbind(dataset7, dataset6) # put them together
names <- row.names(dataset7) #take the names of the genes that we want to keep
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
#expression gene higher than 1000
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
setwd('C:/Users/PAULA CAMARGO ROMERA/OneDrive - Uppsala universitet/Master/Uppsala University/Courses/Year 2/Thesis/processing/cellphonedb_results/normalisation_division')
file <- "matrix_cellphonedb_normalisation_division.txt"
write.table(data, file, append = FALSE, sep = "\t", dec = ".",
            row.names = TRUE, col.names = TRUE)

#metadata normalisation
setwd('C:/Users/PAULA CAMARGO ROMERA/OneDrive - Uppsala universitet/Master/Uppsala University/Courses/Year 2/Thesis/processing/cellphonedb_results/normalisation_division')
x1<-rep(c("AC-like"),each=26101) #glioblastoma AC-like 26101
x2<-rep(c("MES-like"),each=14981) #glioblastoma MES-like 14981
x3<-rep(c("OPC-like"),each=8032) #glioblastoma OPC-like 8032
x4<-rep(c("NPC-like"),each=12294) #glioblastoma NPC-like 12294
x5<-rep(c("microglial cell"),each=26165) #microglial cell (TAM-MG) 26165
x6<-rep(c("macrophage"),each=36356) #macrophages (TAM-BDM) 36356
l1 <- str_split(d1.subset@assays$SCT@data@Dimnames[[2]], '"')
dd1 <- c(l1[1:26101])
dd2 <- c(l1[26102:41082])
dd3 <- c(l1[41083:49114])
dd4 <- c(l1[49115:61408])
dd5 <- c(l1[61409:87573])
dd6 <- c(l1[87574:123929])
df1 <- data.frame(unlist(dd1),x1)
df2 <- data.frame(unlist(dd2),x2)
df3 <- data.frame(unlist(dd3),x3)
df4 <- data.frame(unlist(dd4),x4)
df5 <- data.frame(unlist(dd5),x5)
df6 <- data.frame(unlist(dd6),x6)

names(df1) <- c('Cell', 'cell_type')
names(df2) <- c('Cell', 'cell_type')
names(df3) <- c('Cell', 'cell_type')
names(df4) <- c('Cell', 'cell_type')
names(df5) <- c('Cell', 'cell_type')
names(df6) <- c('Cell', 'cell_type')
df7 <- rbind(df1, df2) # put them together
df7 <- rbind(df7, df3) # put them together
df7 <- rbind(df7, df4) # put them together
df7 <- rbind(df7, df5) # put them together
df7 <- rbind(df7, df6) # put them together
write.table(df7,"meta_cellphone_normalisation_division.txt", row.names = FALSE, sep = "\t")

#post commands
#python commands
#with open('C:/Users/PAULA CAMARGO ROMERA/OneDrive - Uppsala universitet/Master/Uppsala University/Courses/Year 2/Thesis/processing/cellphonedb_results/normalisation_division/meta_cellphone_normalisation_division.txt', 'r', encoding='utf-8') as input_file:
#  lines = input_file.readlines()
#print(lines)
#with open('C:/Users/PAULA CAMARGO ROMERA/OneDrive - Uppsala universitet/Master/Uppsala University/Courses/Year 2/Thesis/processing/cellphonedb_results/normalisation_division/meta_normalisation_division.txt', 'w', encoding='utf-8') as output_file:
#  for line in lines:
#  output_file.write(line.replace('microglial\tcell', 'microglial_cell'))

#linux commands wsl
#sed -e 's/"//g' meta_normalisation_division.txt > meta_normalisation_division2.txt

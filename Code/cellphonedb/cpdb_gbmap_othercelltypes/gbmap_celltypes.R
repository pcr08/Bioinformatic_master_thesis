#load packages
library(Seurat)
library(SingleCellSignalR)
require(Matrix)
require(stringr)

#load matrix
setwd('C:/Users/PAULA CAMARGO ROMERA/OneDrive - Uppsala universitet/Master/Uppsala University/Courses/Year 2/Thesis/processing')
data_newexp <- readRDS('data_newexp.rds')
#four stages (MES/NPC/OPC/AC) of glioblastoma
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


library(SCOPfunctions)
data <- utils_big_as.matrix(mmatrix)
setwd('C:/Users/PAULA CAMARGO ROMERA/OneDrive - Uppsala universitet/Master/Uppsala University/Courses/Year 2/Thesis/processing/cellphonedb_results/gbmap_celltypes')
file <- "matrix_cellphonedb_normalisation_othercelltypes.txt"
write.table(data, file, append = FALSE, sep = "\t", dec = ".",
            row.names = TRUE, col.names = TRUE)

#metadata normalisation
setwd('C:/Users/PAULA CAMARGO ROMERA/OneDrive - Uppsala universitet/Master/Uppsala University/Courses/Year 2/Thesis/processing/cellphonedb_results/gbmap_celltypes')
x1<-rep(c("malignant cell"),each=61408) #malignant cell 61408
x2<-rep(c("mast cell"),each=373) #mast cell 373
x3<-rep(c("endothelial cell"),each=673) #endothelial cell 673
x4<-rep(c("astrocyte"),each=173) #astrocyte 173
x5<-rep(c("oligodendrocyte"),each=12481) #oligodendrocyte 12481
x6<-rep(c("dendritic cell"),each=3961) #dendritic cell 3961
x7<-rep(c("neuron"),each=22) #neuron 22
x8<-rep(c("monocyte"),each=14215) #monocyte 14215
x9<-rep(c("natural killer cell"),each=2489) #natural killer cell 2489
x10<-rep(c("radial glial cell"),each=2807) #radial glial cell 2807
x11<-rep(c("mural cell"),each=1418) #mural cell 1418

l1 <- str_split(se_combined.subset@assays$SCT@data@Dimnames[[2]], '"')
dd1 <- c(l1[1:61408])
dd2 <- c(l1[61409:61781])
dd3 <- c(l1[61782:62454])
dd4 <- c(l1[62455:62627])
dd5 <- c(l1[62628:75108])
dd6 <- c(l1[75109:79069])
dd7 <- c(l1[79070:79091])
dd8 <- c(l1[79092:93306])
dd9 <- c(l1[93307:95795])
dd10 <- c(l1[95796:98602])
dd11 <- c(l1[98603:100020])



df1 <- data.frame(unlist(dd1),x1)
df2 <- data.frame(unlist(dd2),x2)
df3 <- data.frame(unlist(dd3),x3)
df4 <- data.frame(unlist(dd4),x4)
df5 <- data.frame(unlist(dd5),x5)
df6 <- data.frame(unlist(dd6),x6)
df7 <- data.frame(unlist(dd7),x7)
df8 <- data.frame(unlist(dd8),x8)
df9 <- data.frame(unlist(dd9),x9)
df10 <- data.frame(unlist(dd10),x10)
df11 <- data.frame(unlist(dd11),x11)

names(df1) <- c('Cell', 'cell_type')
names(df2) <- c('Cell', 'cell_type')
names(df3) <- c('Cell', 'cell_type')
names(df4) <- c('Cell', 'cell_type')
names(df5) <- c('Cell', 'cell_type')
names(df6) <- c('Cell', 'cell_type')
names(df7) <- c('Cell', 'cell_type')
names(df8) <- c('Cell', 'cell_type')
names(df9) <- c('Cell', 'cell_type')
names(df10) <- c('Cell', 'cell_type')
names(df11) <- c('Cell', 'cell_type')

df12 <- rbind(df1, df2) # put them together
df12 <- rbind(df12, df3) # put them together
df12 <- rbind(df12, df4) # put them together
df12 <- rbind(df12, df5) # put them together
df12 <- rbind(df12, df6) # put them together
df12 <- rbind(df12, df7) # put them together
df12 <- rbind(df12, df8) # put them together
df12 <- rbind(df12, df9) # put them together
df12 <- rbind(df12, df10) # put them together
df12 <- rbind(df12, df11) # put them together
write.table(df12,"meta_cellphone_normalisation_othercelltypes.txt", row.names = FALSE, sep = "\t")



#post commands
#python commands
#with open('C:/Users/PAULA CAMARGO ROMERA/OneDrive - Uppsala universitet/Master/Uppsala University/Courses/Year 2/Thesis/processing/cellphonedb_results/normalisation_division/meta_cellphone_normalisation_division.txt', 'r', encoding='utf-8') as input_file:
#  lines = input_file.readlines()
#print(lines)
#with open('C:/Users/PAULA CAMARGO ROMERA/OneDrive - Uppsala universitet/Master/Uppsala University/Courses/Year 2/Thesis/processing/cellphonedb_results/normalisation_division/meta_normalisation_division.txt', 'w', encoding='utf-8') as output_file:
#  for line in lines:
#  output_file.write(line.replace('malignant\tcell', 'malignant_cell'))
#  output_file.write(line.replace('mast\tcell', 'mast_cell'))
#  output_file.write(line.replace('endothelial\tcell', 'endothelial_cell'))
#  output_file.write(line.replace('dendritic\tcell', 'dendritic_cell'))
#  output_file.write(line.replace('natural\tkiller\tcell', 'natural_killer_cell'))
#  output_file.write(line.replace('radial\tglial\tcell', 'radial_glial_cell'))
#  output_file.write(line.replace('mural\tcell', 'mural_cell'))

#linux commands wsl
#sed -e 's/"//g' meta_normalisation_division.txt > meta_normalisation_division2.txt

library("readxl")
d1 <- read_excel("C:/Users/PAULA CAMARGO ROMERA/OneDrive - Uppsala universitet/Master/Uppsala University/Courses/Year 2/Thesis/processing/cellphonedb_results/gbmap_celltypes/out/significant_means2_radial_malig.xlsx")
l1 <- d1[[2]]

v1 <- ""
for (i in l1){
print(i)
  v1 <- paste(v1, i, sep = '","')
}
v1
  
  
  
  
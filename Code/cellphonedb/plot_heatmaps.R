library(pheatmap)
heatmaps_plot = function(meta_file, pvalues_file, count_filename, log_filename, count_network_filename, interaction_count_filename, count_network_separator, interaction_count_separator, show_rownames = T, show_colnames = T,
                         scale="none", cluster_cols = T,border_color='white', cluster_rows = T, fontsize_row=11,
                         fontsize_col = 11, main = '',treeheight_row=0, family='Arial', treeheight_col = 0,
                         col1 = "dodgerblue4", col2 = 'peachpuff', col3 = 'deeppink4', meta_sep='\t', pvalues_sep='\t', pvalue=0.05){
  #######   Network
  
  meta = read.csv(meta_file, comment.char = '', sep=meta_sep)
  
  all_intr = read.table(pvalues_file, header=T, stringsAsFactors = F, sep=pvalues_sep, comment.char = '', check.names = F)
  intr_pairs = all_intr$interacting_pair
  all_intr = all_intr[,-c(1:11)]
  
  
  split_sep = '\\|'
  join_sep = '|'
  
  pairs1_all = unique(meta[,2])
  
  pairs1 = c()
  for (i in 1:length(pairs1_all))
    for (j in 1:length(pairs1_all))
      pairs1 = c(pairs1,paste(pairs1_all[i],pairs1_all[j],sep=join_sep))
  
  all_count = matrix(ncol=3)
  colnames(all_count) = c('SOURCE','TARGET','count')
  count1 = c()
  for(i in 1:length(pairs1))
  {
    p1 = strsplit(pairs1[i], split_sep)[[1]][1]
    p2 = strsplit(pairs1[i], split_sep)[[1]][2]
    
    n1 = intr_pairs[which(all_intr[,pairs1[i]]<=pvalue)]
    pairs_rev = paste(p2, p1, sep=join_sep)
    n2 = intr_pairs[which(all_intr[,pairs_rev]<=pvalue)]
    
    if(p1!=p2)
      count1 = length(unique(n1))+length(unique(n2))
    else
      count1 = length(unique(n1))
    
    new_count = c(p1,p2,count1)
    names(new_count) = c('SOURCE','TARGET','count')
    all_count = rbind(all_count, new_count)
  }
  
  all_count = all_count[-1,]
  write.table(all_count, count_network_filename, sep=count_network_separator, quote=F, row.names = F)
  
  #######   count interactions
  
  count1 = c()
  for(i in 1:length(pairs1))
  {
    p1 = strsplit(pairs1[i], split_sep)[[1]][1]
    p2 = strsplit(pairs1[i], split_sep)[[1]][2]
    
    n1 = intr_pairs[which(all_intr[,pairs1[i]]<=pvalue)]
    
    pairs_rev = paste(p2, p1, sep=join_sep)
    n2 = intr_pairs[which(all_intr[,pairs_rev]<=pvalue)]
    if(p1!=p2)
      count1 = c(count1,length(unique(n1))+length(unique(n2)))
    else
      count1 = c(count1,length(unique(n1)))
    
  }
  
  if (any(count1)>0)
  {
    count_matrix = matrix(count1, nrow=length(unique(meta[,2])), ncol=length(unique(meta[,2])))
    rownames(count_matrix)= unique(meta[,2])
    colnames(count_matrix)= unique(meta[,2])
    
    all_sum = rowSums(count_matrix)
    all_sum = cbind(names(all_sum), all_sum)
    write.table(all_sum, file=interaction_count_filename, quote=F, sep=count_network_separator, row.names=F)
    
    col.heatmap <- colorRampPalette(c(col1,col2,col3 ))( 1000 )
    
    pheatmap(count_matrix, show_rownames = show_rownames, show_colnames = show_colnames, scale=scale, cluster_cols = cluster_cols,
             border_color=border_color, cluster_rows = cluster_rows, fontsize_row = fontsize_row, fontsize_col = fontsize_col,
             main = main, treeheight_row = treeheight_row, family = family,color = col.heatmap, treeheight_col = treeheight_col, filename = count_filename)
    
    pheatmap(log(count_matrix+1), show_rownames = show_rownames, show_colnames = show_colnames, scale=scale, cluster_cols = cluster_cols,
             border_color=border_color, cluster_rows = cluster_rows, fontsize_row = fontsize_row, fontsize_col = fontsize_col,
             main = main, treeheight_row = treeheight_row, family = family,color = col.heatmap, treeheight_col = treeheight_col, filename = log_filename)
  } else {
    stop("There are no significant results using p-value of: ", pvalue, call.=FALSE)
  }
}

#gmap microglia-glioblastoma
#setwd('C:/Users/PAULA CAMARGO ROMERA/OneDrive - Uppsala universitet/Master/Uppsala University/Courses/Year 2/Thesis/processing/cellphonedb_results/normalisation')
#heatmaps_plot(meta_file = "meta_data.txt", pvalues_file = "pvalues.txt", count_network_filename = "network.txt", count_network_separator = "\t", interaction_count_filename = "interaction_count.txt", count_filename = "heatmap_count.pdf", log_filename = "heatmap_log_count.pdf")

#division
#setwd('C:/Users/PAULA CAMARGO ROMERA/OneDrive - Uppsala universitet/Master/Uppsala University/Courses/Year 2/Thesis/processing/cellphonedb_results/normalisation_division')
#heatmaps_plot(meta_file = "meta_normalisation_division.txt", pvalues_file = "pvalues.txt", count_network_filename = "network.txt", count_network_separator = "\t", interaction_count_filename = "interaction_count.txt", count_filename = "heatmap_count.pdf", log_filename = "heatmap_log_count.pdf")

#our data & gbmap
#setwd('C:/Users/PAULA CAMARGO ROMERA/OneDrive - Uppsala universitet/Master/Uppsala University/Courses/Year 2/Thesis/processing/our_data/cellphonedb_results')
#heatmaps_plot(meta_file = "meta_normalisation_ourdata2.txt", pvalues_file = "pvalues.txt", count_network_filename = "network.txt", count_network_separator = "\t", interaction_count_filename = "interaction_count.txt", count_filename = "heatmap_count.pdf", log_filename = "heatmap_log_count.pdf")


#our data patients
#U3005
setwd('D:/new_boot/neoplastic_U3005/out')
heatmaps_plot(meta_file = "../meta_cellp_normalisation_ourdata_bootstrap_meta_cellphone_normalisation_ourdata_bootstrap_indneoplastic_U3005_2.txt", pvalues_file = "pvalues.txt", count_network_filename = "network.txt", count_network_separator = "\t", interaction_count_filename = "interaction_count.txt", count_filename = "heatmap_count.pdf", log_filename = "heatmap_log_count.pdf")

#U3008
setwd('D:/new_boot/neoplastic_U3008/out')
heatmaps_plot(meta_file = "../meta_cellp_normalisation_ourdata_bootstrap_meta_cellphone_normalisation_ourdata_bootstrap_indneoplastic_U3008_2.txt", pvalues_file = "pvalues.txt", count_network_filename = "network.txt", count_network_separator = "\t", interaction_count_filename = "interaction_count.txt", count_filename = "heatmap_count.pdf", log_filename = "heatmap_log_count.pdf")


#U3013
setwd('D:/new_boot/neoplastic_U3013/out')
heatmaps_plot(meta_file = "../meta_cellp_normalisation_ourdata_bootstrap_meta_cellphone_normalisation_ourdata_bootstrap_indneoplastic_U3013_2.txt", pvalues_file = "pvalues.txt", count_network_filename = "network.txt", count_network_separator = "\t", interaction_count_filename = "interaction_count.txt", count_filename = "heatmap_count.pdf", log_filename = "heatmap_log_count.pdf")


#U3016
setwd('D:/new_boot/neoplastic_U3016/out')
heatmaps_plot(meta_file = "../meta_cellp_normalisation_ourdata_bootstrap_meta_cellphone_normalisation_ourdata_bootstrap_indneoplastic_U3016_2.txt", pvalues_file = "pvalues.txt", count_network_filename = "network.txt", count_network_separator = "\t", interaction_count_filename = "interaction_count.txt", count_filename = "heatmap_count.pdf", log_filename = "heatmap_log_count.pdf")


#U3017
setwd('D:/new_boot/neoplastic_U3017/out')
heatmaps_plot(meta_file = "../meta_cellp_normalisation_ourdata_bootstrap_meta_cellphone_normalisation_ourdata_bootstrap_indneoplastic_U3017_2.txt", pvalues_file = "pvalues.txt", count_network_filename = "network.txt", count_network_separator = "\t", interaction_count_filename = "interaction_count.txt", count_filename = "heatmap_count.pdf", log_filename = "heatmap_log_count.pdf")

#U3020
setwd('D:/new_boot/neoplastic_U3020/out')
heatmaps_plot(meta_file = "../meta_cellp_normalisation_ourdata_bootstrap_meta_cellphone_normalisation_ourdata_bootstrap_indneoplastic_U3020_2.txt", pvalues_file = "pvalues.txt", count_network_filename = "network.txt", count_network_separator = "\t", interaction_count_filename = "interaction_count.txt", count_filename = "heatmap_count.pdf", log_filename = "heatmap_log_count.pdf")


#U3021
setwd('D:/new_boot/neoplastic_U3021/out')
heatmaps_plot(meta_file = "../meta_cellp_normalisation_ourdata_bootstrap_meta_cellphone_normalisation_ourdata_bootstrap_indneoplastic_U3021_2.txt", pvalues_file = "pvalues.txt", count_network_filename = "network.txt", count_network_separator = "\t", interaction_count_filename = "interaction_count.txt", count_filename = "heatmap_count.pdf", log_filename = "heatmap_log_count.pdf")


#U3024
setwd('D:/new_boot/neoplastic_U3024/out')
heatmaps_plot(meta_file = "../meta_cellp_normalisation_ourdata_bootstrap_meta_cellphone_normalisation_ourdata_bootstrap_indneoplastic_U3024_2.txt", pvalues_file = "pvalues.txt", count_network_filename = "network.txt", count_network_separator = "\t", interaction_count_filename = "interaction_count.txt", count_filename = "heatmap_count.pdf", log_filename = "heatmap_log_count.pdf")


#U3027
setwd('D:/new_boot/neoplastic_U3027/out')
heatmaps_plot(meta_file = "../meta_cellp_normalisation_ourdata_bootstrap_meta_cellphone_normalisation_ourdata_bootstrap_indneoplastic_U3027_2.txt", pvalues_file = "pvalues.txt", count_network_filename = "network.txt", count_network_separator = "\t", interaction_count_filename = "interaction_count.txt", count_filename = "heatmap_count.pdf", log_filename = "heatmap_log_count.pdf")


#U3028
setwd('D:/new_boot/neoplastic_U3028/out')
heatmaps_plot(meta_file = "../meta_cellp_normalisation_ourdata_bootstrap_meta_cellphone_normalisation_ourdata_bootstrap_indneoplastic_U3028_2.txt", pvalues_file = "pvalues.txt", count_network_filename = "network.txt", count_network_separator = "\t", interaction_count_filename = "interaction_count.txt", count_filename = "heatmap_count.pdf", log_filename = "heatmap_log_count.pdf")


#U3031
setwd('D:/new_boot/neoplastic_U3031/out')
heatmaps_plot(meta_file = "../meta_cellp_normalisation_ourdata_bootstrap_meta_cellphone_normalisation_ourdata_bootstrap_indneoplastic_U3031_2.txt", pvalues_file = "pvalues.txt", count_network_filename = "network.txt", count_network_separator = "\t", interaction_count_filename = "interaction_count.txt", count_filename = "heatmap_count.pdf", log_filename = "heatmap_log_count.pdf")


#U3033
setwd('D:/new_boot/neoplastic_U3033/out')
heatmaps_plot(meta_file = "../meta_cellp_normalisation_ourdata_bootstrap_meta_cellphone_normalisation_ourdata_bootstrap_indneoplastic_U3033_2.txt", pvalues_file = "pvalues.txt", count_network_filename = "network.txt", count_network_separator = "\t", interaction_count_filename = "interaction_count.txt", count_filename = "heatmap_count.pdf", log_filename = "heatmap_log_count.pdf")


#U3035
setwd('D:/new_boot/neoplastic_U3035/out')
heatmaps_plot(meta_file = "../meta_cellp_normalisation_ourdata_bootstrap_meta_cellphone_normalisation_ourdata_bootstrap_indneoplastic_U3035_2.txt", pvalues_file = "pvalues.txt", count_network_filename = "network.txt", count_network_separator = "\t", interaction_count_filename = "interaction_count.txt", count_filename = "heatmap_count.pdf", log_filename = "heatmap_log_count.pdf")


#U3039
setwd('D:/new_boot/neoplastic_U3039/out')
heatmaps_plot(meta_file = "../meta_cellp_normalisation_ourdata_bootstrap_meta_cellphone_normalisation_ourdata_bootstrap_indneoplastic_U3039_2.txt", pvalues_file = "pvalues.txt", count_network_filename = "network.txt", count_network_separator = "\t", interaction_count_filename = "interaction_count.txt", count_filename = "heatmap_count.pdf", log_filename = "heatmap_log_count.pdf")


#U3047
setwd('D:/new_boot/neoplastic_U3047/out')
heatmaps_plot(meta_file = "../meta_cellp_normalisation_ourdata_bootstrap_meta_cellphone_normalisation_ourdata_bootstrap_indneoplastic_U3047_2.txt", pvalues_file = "pvalues.txt", count_network_filename = "network.txt", count_network_separator = "\t", interaction_count_filename = "interaction_count.txt", count_filename = "heatmap_count.pdf", log_filename = "heatmap_log_count.pdf")


#U3051
setwd('D:/new_boot/neoplastic_U3051/out')
heatmaps_plot(meta_file = "../meta_cellp_normalisation_ourdata_bootstrap_meta_cellphone_normalisation_ourdata_bootstrap_indneoplastic_U3051_2.txt", pvalues_file = "pvalues.txt", count_network_filename = "network.txt", count_network_separator = "\t", interaction_count_filename = "interaction_count.txt", count_filename = "heatmap_count.pdf", log_filename = "heatmap_log_count.pdf")

#U3053
setwd('D:/new_boot/neoplastic_U3053/out')
heatmaps_plot(meta_file = "../meta_cellp_normalisation_ourdata_bootstrap_meta_cellphone_normalisation_ourdata_bootstrap_indneoplastic_U3053_2.txt", pvalues_file = "pvalues.txt", count_network_filename = "network.txt", count_network_separator = "\t", interaction_count_filename = "interaction_count.txt", count_filename = "heatmap_count.pdf", log_filename = "heatmap_log_count.pdf")


#U3054
setwd('D:/new_boot/neoplastic_U3054/out')
heatmaps_plot(meta_file = "../meta_cellp_normalisation_ourdata_bootstrap_meta_cellphone_normalisation_ourdata_bootstrap_indneoplastic_U3054_2.txt", pvalues_file = "pvalues.txt", count_network_filename = "network.txt", count_network_separator = "\t", interaction_count_filename = "interaction_count.txt", count_filename = "heatmap_count.pdf", log_filename = "heatmap_log_count.pdf")


#U3065
setwd('D:/new_boot/neoplastic_U3065/out')
heatmaps_plot(meta_file = "../meta_cellp_normalisation_ourdata_bootstrap_meta_cellphone_normalisation_ourdata_bootstrap_indneoplastic_U3065_2.txt", pvalues_file = "pvalues.txt", count_network_filename = "network.txt", count_network_separator = "\t", interaction_count_filename = "interaction_count.txt", count_filename = "heatmap_count.pdf", log_filename = "heatmap_log_count.pdf")


#U3071
setwd('D:/new_boot/neoplastic_U3071/out')
heatmaps_plot(meta_file = "../meta_cellp_normalisation_ourdata_bootstrap_meta_cellphone_normalisation_ourdata_bootstrap_indneoplastic_U3071_2.txt", pvalues_file = "pvalues.txt", count_network_filename = "network.txt", count_network_separator = "\t", interaction_count_filename = "interaction_count.txt", count_filename = "heatmap_count.pdf", log_filename = "heatmap_log_count.pdf")


#U3078
setwd('D:/new_boot/neoplastic_U3078/out')
heatmaps_plot(meta_file = "../meta_cellp_normalisation_ourdata_bootstrap_meta_cellphone_normalisation_ourdata_bootstrap_indneoplastic_U3078_2.txt", pvalues_file = "pvalues.txt", count_network_filename = "network.txt", count_network_separator = "\t", interaction_count_filename = "interaction_count.txt", count_filename = "heatmap_count.pdf", log_filename = "heatmap_log_count.pdf")

#U3082
setwd('D:/new_boot/neoplastic_U3082/out')
heatmaps_plot(meta_file = "../meta_cellp_normalisation_ourdata_bootstrap_meta_cellphone_normalisation_ourdata_bootstrap_indneoplastic_U3082_2.txt", pvalues_file = "pvalues.txt", count_network_filename = "network.txt", count_network_separator = "\t", interaction_count_filename = "interaction_count.txt", count_filename = "heatmap_count.pdf", log_filename = "heatmap_log_count.pdf")


#U3084
setwd('D:/new_boot/neoplastic_U3084/out')
heatmaps_plot(meta_file = "../meta_cellp_normalisation_ourdata_bootstrap_meta_cellphone_normalisation_ourdata_bootstrap_indneoplastic_U3084_2.txt", pvalues_file = "pvalues.txt", count_network_filename = "network.txt", count_network_separator = "\t", interaction_count_filename = "interaction_count.txt", count_filename = "heatmap_count.pdf", log_filename = "heatmap_log_count.pdf")


#U3086
setwd('D:/new_boot/neoplastic_U3086/out')
heatmaps_plot(meta_file = "../meta_cellp_normalisation_ourdata_bootstrap_meta_cellphone_normalisation_ourdata_bootstrap_indneoplastic_U3086_2.txt", pvalues_file = "pvalues.txt", count_network_filename = "network.txt", count_network_separator = "\t", interaction_count_filename = "interaction_count.txt", count_filename = "heatmap_count.pdf", log_filename = "heatmap_log_count.pdf")


#U3088
setwd('D:/new_boot/neoplastic_U3088/out')
heatmaps_plot(meta_file = "../meta_cellp_normalisation_ourdata_bootstrap_meta_cellphone_normalisation_ourdata_bootstrap_indneoplastic_U3088_2.txt", pvalues_file = "pvalues.txt", count_network_filename = "network.txt", count_network_separator = "\t", interaction_count_filename = "interaction_count.txt", count_filename = "heatmap_count.pdf", log_filename = "heatmap_log_count.pdf")


#U3093
setwd('D:/new_boot/neoplastic_U3093/out')
heatmaps_plot(meta_file = "../meta_cellp_normalisation_ourdata_bootstrap_meta_cellphone_normalisation_ourdata_bootstrap_indneoplastic_U3093_2.txt", pvalues_file = "pvalues.txt", count_network_filename = "network.txt", count_network_separator = "\t", interaction_count_filename = "interaction_count.txt", count_filename = "heatmap_count.pdf", log_filename = "heatmap_log_count.pdf")


#U3101
setwd('D:/new_boot/neoplastic_U3101/out')
heatmaps_plot(meta_file = "../meta_cellp_normalisation_ourdata_bootstrap_meta_cellphone_normalisation_ourdata_bootstrap_indneoplastic_U3101_2.txt", pvalues_file = "pvalues.txt", count_network_filename = "network.txt", count_network_separator = "\t", interaction_count_filename = "interaction_count.txt", count_filename = "heatmap_count.pdf", log_filename = "heatmap_log_count.pdf")


#U3110
setwd('D:/new_boot/neoplastic_U3110/out')
heatmaps_plot(meta_file = "../meta_cellp_normalisation_ourdata_bootstrap_meta_cellphone_normalisation_ourdata_bootstrap_indneoplastic_U3110_2.txt", pvalues_file = "pvalues.txt", count_network_filename = "network.txt", count_network_separator = "\t", interaction_count_filename = "interaction_count.txt", count_filename = "heatmap_count.pdf", log_filename = "heatmap_log_count.pdf")


#U3117
setwd('D:/new_boot/neoplastic_U3117/out')
heatmaps_plot(meta_file = "../meta_cellp_normalisation_ourdata_bootstrap_meta_cellphone_normalisation_ourdata_bootstrap_indneoplastic_U3117_2.txt", pvalues_file = "pvalues.txt", count_network_filename = "network.txt", count_network_separator = "\t", interaction_count_filename = "interaction_count.txt", count_filename = "heatmap_count.pdf", log_filename = "heatmap_log_count.pdf")

#U3118
setwd('D:/Results/cellphonedb/cpdb_new_boot_ind/neoplastic_U3118/out')
heatmaps_plot(meta_file = "../meta_cellp_normalisation_ourdata_bootstrap_meta_cellphone_normalisation_ourdata_bootstrap_indneoplastic_U3118_2.txt", pvalues_file = "pvalues.txt", count_network_filename = "network.txt", count_network_separator = "\t", interaction_count_filename = "interaction_count.txt", count_filename = "heatmap_count.pdf", log_filename = "heatmap_log_count.pdf")


#U3123
setwd('D:/Results/cellphonedb/cpdb_new_boot_ind/neoplastic_U3123/out')
heatmaps_plot(meta_file = "../meta_cellp_normalisation_ourdata_bootstrap_meta_cellphone_normalisation_ourdata_bootstrap_indneoplastic_U3123_2.txt", pvalues_file = "pvalues.txt", count_network_filename = "network.txt", count_network_separator = "\t", interaction_count_filename = "interaction_count.txt", count_filename = "heatmap_count.pdf", log_filename = "heatmap_log_count.pdf")

#U3137
setwd('D:/Results/cellphonedb/cpdb_new_boot_ind/neoplastic_U3137/out')
heatmaps_plot(meta_file = "../meta_cellp_normalisation_ourdata_bootstrap_meta_cellphone_normalisation_ourdata_bootstrap_indneoplastic_U3137_2.txt", pvalues_file = "pvalues.txt", count_network_filename = "network.txt", count_network_separator = "\t", interaction_count_filename = "interaction_count.txt", count_filename = "heatmap_count.pdf", log_filename = "heatmap_log_count.pdf")

#U3167
setwd('D:/Results/cellphonedb/cpdb_new_boot_ind/neoplastic_U3167/out')
heatmaps_plot(meta_file = "../meta_cellp_normalisation_ourdata_bootstrap_meta_cellphone_normalisation_ourdata_bootstrap_indneoplastic_U3167_2.txt", pvalues_file = "pvalues.txt", count_network_filename = "network.txt", count_network_separator = "\t", interaction_count_filename = "interaction_count.txt", count_filename = "heatmap_count.pdf", log_filename = "heatmap_log_count.pdf")


#U3169
setwd('D:/Results/cellphonedb/cpdb_new_boot_ind/neoplastic_U3169/out')
heatmaps_plot(meta_file = "../meta_cellp_normalisation_ourdata_bootstrap_meta_cellphone_normalisation_ourdata_bootstrap_indneoplastic_U3169_2.txt", pvalues_file = "pvalues.txt", count_network_filename = "network.txt", count_network_separator = "\t", interaction_count_filename = "interaction_count.txt", count_filename = "heatmap_count.pdf", log_filename = "heatmap_log_count.pdf")

#U3179
setwd('D:/Results/cellphonedb/cpdb_new_boot_ind/neoplastic_U3179/out')
heatmaps_plot(meta_file = "../meta_cellp_normalisation_ourdata_bootstrap_meta_cellphone_normalisation_ourdata_bootstrap_indneoplastic_U3179_2.txt", pvalues_file = "pvalues.txt", count_network_filename = "network.txt", count_network_separator = "\t", interaction_count_filename = "interaction_count.txt", count_filename = "heatmap_count.pdf", log_filename = "heatmap_log_count.pdf")


#U3180
setwd('D:/new_boot/neoplastic_U3180/out')
heatmaps_plot(meta_file = "../meta_cellp_normalisation_ourdata_bootstrap_meta_cellphone_normalisation_ourdata_bootstrap_indneoplastic_U3180_2.txt", pvalues_file = "pvalues.txt", count_network_filename = "network.txt", count_network_separator = "\t", interaction_count_filename = "interaction_count.txt", count_filename = "heatmap_count.pdf", log_filename = "heatmap_log_count.pdf")


#U3182
setwd('D:/new_boot/neoplastic_U3182/out')
heatmaps_plot(meta_file = "../meta_cellp_normalisation_ourdata_bootstrap_meta_cellphone_normalisation_ourdata_bootstrap_indneoplastic_U3180_2.txt", pvalues_file = "pvalues.txt", count_network_filename = "network.txt", count_network_separator = "\t", interaction_count_filename = "interaction_count.txt", count_filename = "heatmap_count.pdf", log_filename = "heatmap_log_count.pdf")


#U3202
setwd('D:/Results/cellphonedb/cpdb_new_boot_ind/neoplastic_U3202/out')
heatmaps_plot(meta_file = "../meta_cellp_normalisation_ourdata_bootstrap_meta_cellphone_normalisation_ourdata_bootstrap_indneoplastic_U3202_2.txt", pvalues_file = "pvalues.txt", count_network_filename = "network.txt", count_network_separator = "\t", interaction_count_filename = "interaction_count.txt", count_filename = "heatmap_count.pdf", log_filename = "heatmap_log_count.pdf")


#U3213
setwd('D:/Results/cellphonedb/cpdb_new_boot_ind/neoplastic_U3213/out')
heatmaps_plot(meta_file = "../meta_cellp_normalisation_ourdata_bootstrap_meta_cellphone_normalisation_ourdata_bootstrap_indneoplastic_U3213_2.txt", pvalues_file = "pvalues.txt", count_network_filename = "network.txt", count_network_separator = "\t", interaction_count_filename = "interaction_count.txt", count_filename = "heatmap_count.pdf", log_filename = "heatmap_log_count.pdf")


#U3220
setwd('D:/Results/cellphonedb/cpdb_new_boot_ind/neoplastic_U3220/out')
heatmaps_plot(meta_file = "../meta_cellp_normalisation_ourdata_bootstrap_meta_cellphone_normalisation_ourdata_bootstrap_indneoplastic_U3220_2.txt", pvalues_file = "pvalues.txt", count_network_filename = "network.txt", count_network_separator = "\t", interaction_count_filename = "interaction_count.txt", count_filename = "heatmap_count.pdf", log_filename = "heatmap_log_count.pdf")


#U3226
setwd('D:/Results/cellphonedb/cpdb_new_boot_ind/neoplastic_U3226/out')
heatmaps_plot(meta_file = "../meta_cellp_normalisation_ourdata_bootstrap_meta_cellphone_normalisation_ourdata_bootstrap_indneoplastic_U3226_2.txt", pvalues_file = "pvalues.txt", count_network_filename = "network.txt", count_network_separator = "\t", interaction_count_filename = "interaction_count.txt", count_filename = "heatmap_count.pdf", log_filename = "heatmap_log_count.pdf")


#U3230
setwd('D:/Results/cellphonedb/cpdb_new_boot_ind/neoplastic_U3230/out')
heatmaps_plot(meta_file = "../meta_cellp_normalisation_ourdata_bootstrap_meta_cellphone_normalisation_ourdata_bootstrap_indneoplastic_U3230_2.txt", pvalues_file = "pvalues.txt", count_network_filename = "network.txt", count_network_separator = "\t", interaction_count_filename = "interaction_count.txt", count_filename = "heatmap_count.pdf", log_filename = "heatmap_log_count.pdf")


#U3233
setwd('D:/Results/cellphonedb/cpdb_new_boot_ind/neoplastic_U3233/out')
heatmaps_plot(meta_file = "../meta_cellp_normalisation_ourdata_bootstrap_meta_cellphone_normalisation_ourdata_bootstrap_indneoplastic_U3233_2.txt", pvalues_file = "pvalues.txt", count_network_filename = "network.txt", count_network_separator = "\t", interaction_count_filename = "interaction_count.txt", count_filename = "heatmap_count.pdf", log_filename = "heatmap_log_count.pdf")


#U3257
setwd('D:/Results/cellphonedb/cpdb_new_boot_ind/neoplastic_U3257/out')
heatmaps_plot(meta_file = "../meta_cellp_normalisation_ourdata_bootstrap_meta_cellphone_normalisation_ourdata_bootstrap_indneoplastic_U3257_2.txt", pvalues_file = "pvalues.txt", count_network_filename = "network.txt", count_network_separator = "\t", interaction_count_filename = "interaction_count.txt", count_filename = "heatmap_count.pdf", log_filename = "heatmap_log_count.pdf")


#U3275
setwd('D:/Results/cellphonedb/cpdb_new_boot_ind/neoplastic_U3275/out')
heatmaps_plot(meta_file = "../meta_cellp_normalisation_ourdata_bootstrap_meta_cellphone_normalisation_ourdata_bootstrap_indneoplastic_U3275_2.txt", pvalues_file = "pvalues.txt", count_network_filename = "network.txt", count_network_separator = "\t", interaction_count_filename = "interaction_count.txt", count_filename = "heatmap_count.pdf", log_filename = "heatmap_log_count.pdf")


#U3279
setwd('D:/Results/cellphonedb/cpdb_new_boot_ind/neoplastic_U3279/out')
heatmaps_plot(meta_file = "../meta_cellp_normalisation_ourdata_bootstrap_meta_cellphone_normalisation_ourdata_bootstrap_indneoplastic_U3279_2.txt", pvalues_file = "pvalues.txt", count_network_filename = "network.txt", count_network_separator = "\t", interaction_count_filename = "interaction_count.txt", count_filename = "heatmap_count.pdf", log_filename = "heatmap_log_count.pdf")


#U3281
setwd('D:/Results/cellphonedb/cpdb_new_boot_ind/neoplastic_U3281/out')
heatmaps_plot(meta_file = "../meta_cellp_normalisation_ourdata_bootstrap_meta_cellphone_normalisation_ourdata_bootstrap_indneoplastic_U3281_2.txt", pvalues_file = "pvalues.txt", count_network_filename = "network.txt", count_network_separator = "\t", interaction_count_filename = "interaction_count.txt", count_filename = "heatmap_count.pdf", log_filename = "heatmap_log_count.pdf")


#U3289
setwd('D:/Results/cellphonedb/cpdb_new_boot_ind/neoplastic_U3289/out')
heatmaps_plot(meta_file = "../meta_cellp_normalisation_ourdata_bootstrap_meta_cellphone_normalisation_ourdata_bootstrap_indneoplastic_U3289_2.txt", pvalues_file = "pvalues.txt", count_network_filename = "network.txt", count_network_separator = "\t", interaction_count_filename = "interaction_count.txt", count_filename = "heatmap_count.pdf", log_filename = "heatmap_log_count.pdf")

#U3291
setwd('D:/Results/cellphonedb/cpdb_new_boot_ind/neoplastic_U3291/out')
heatmaps_plot(meta_file = "../meta_cellp_normalisation_ourdata_bootstrap_meta_cellphone_normalisation_ourdata_bootstrap_indneoplastic_U3291_2.txt", pvalues_file = "pvalues.txt", count_network_filename = "network.txt", count_network_separator = "\t", interaction_count_filename = "interaction_count.txt", count_filename = "heatmap_count.pdf", log_filename = "heatmap_log_count.pdf")


#U3299
setwd('D:/Results/cellphonedb/cpdb_new_boot_ind/neoplastic_U3299/out')
heatmaps_plot(meta_file = "../meta_cellp_normalisation_ourdata_bootstrap_meta_cellphone_normalisation_ourdata_bootstrap_indneoplastic_U3299_2.txt", pvalues_file = "pvalues.txt", count_network_filename = "network.txt", count_network_separator = "\t", interaction_count_filename = "interaction_count.txt", count_filename = "heatmap_count.pdf", log_filename = "heatmap_log_count.pdf")


#GMAP glioblastoma vs other cell types
setwd('C:/Users/PAULA CAMARGO ROMERA/OneDrive - Uppsala universitet/Master/Uppsala University/Courses/Year 2/Thesis/processing/cellphonedb_results/gbmap_celltypes/out')
heatmaps_plot(meta_file = "../meta_cellphone_normalisation_othercelltypes2.txt", pvalues_file = "pvalues.txt", count_network_filename = "network.txt", count_network_separator = "\t", interaction_count_filename = "interaction_count.txt", count_filename = "heatmap_count.pdf", log_filename = "heatmap_log_count.pdf")


import os

#metadata transform
directory = 'C:/Users/PAULA CAMARGO ROMERA/OneDrive - Uppsala universitet/Master/Uppsala University/Courses/Year 2/Thesis/processing/cellphonedb_results/gbmap_celltypes/mural_cells/'
print(directory)
for filename in os.listdir(directory):
    #print(filename)
    path = directory + filename
    with open(path, 'r', encoding='utf-8') as input_file:
        lines = input_file.readlines()
        print(lines)
        file = 'C:/Users/PAULA CAMARGO ROMERA/OneDrive - Uppsala universitet/Master/Uppsala University/Courses/Year 2/Thesis/processing/cellphonedb_results/gbmap_celltypes/mural_cells/meta_cellp_othercelltypes_' + filename
        print(file)
        with open(file, 'w', encoding='utf-8') as output_file:
            for line in lines:
                #output_file.write(line.replace('microglial cell', 'microglial_cell'))
                output_file.write(line.replace('"', ''))


#python commands
# with open('C:/Users/PAULA CAMARGO ROMERA/OneDrive - Uppsala universitet/Master/Uppsala University/Courses/Year 2/Thesis/processing/cellphonedb_results/gbmap_celltypes/meta_cellphone_normalisation_ourdata_bootstrap_indneoplastic_U3005.txt', 'r', encoding='utf-8') as input_file:
#     lines = input_file.readlines()
#     print(lines)
#     with open('C:/Users/PAULA CAMARGO ROMERA/OneDrive - Uppsala universitet/Master/Uppsala University/Courses/Year 2/Thesis/processing/cellphonedb_results/gbmap_celltypes/meta_cellphone_normalisation_ourdata_bootstrap_indneoplastic_U3005_2.txt', 'w', encoding='utf-8') as output_file:
#         for line in lines:
#             #output_file.write(line.replace('malignant cell', 'malignant_cell'))
#             #output_file.write(line.replace('mast cell', 'mast_cell'))
#             #output_file.write(line.replace('endothelial cell', 'endothelial_cell'))
#             #output_file.write(line.replace('dendritic cell', 'dendritic_cell'))
#             #output_file.write(line.replace('natural killer cell', 'natural_killer_cell'))
#             #output_file.write(line.replace('radial glial cell', 'radial_glial_cell'))
#             #output_file.write(line.replace('mural cell', 'mural_cell'))

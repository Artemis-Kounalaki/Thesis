import ast
import pandas as pd
import os

# Open the file with proteins having overlapped ids.

os.chdir(os.path.expanduser('~/conserved_gene_order/human_reference'))
with open("overlapped_ids.txt") as file:
    gene_names = ast.literal_eval(file.read())


# Open human blast results with pandas and delete those rows
# that contain the overlapped proteins.

data = pd.read_csv('results_human.txt', sep="\t", header=None)
data=data[~data[0].isin(gene_names)]
data=data[~data[1].isin(gene_names)]
data.to_csv(r'~/conserved_gene_order/human_reference/clean_blast_results_human.txt', header=None, index=None, sep=' ', mode='a')


# HOW MANY UNIQUE IDS NOW REMAIN.

#column_values = data[[0, 1]].values.ravel()
#unique_values =  pd.unique(column_values)
#print(len(unique_values))

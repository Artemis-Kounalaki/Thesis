import ast
import pandas as pd
import os

# Open the file with human proteins having overlapped ids- isoforms-.

os.chdir(os.path.expanduser('~/conserved_gene_order_demo/human_reference'))

with open("overlapped_ids.txt") as file:
    gene_names = ast.literal_eval(file.read())

os.chdir(os.path.expanduser('~/conserved_gene_order_demo/macaca_reference'))

with open("overlapped_ids_macaca.txt") as file:
    gene_names_mac = ast.literal_eval(file.read())

# Open human-macaca blast results with pandas and delete those rows
# that contain the overlapped proteins.

data = pd.read_csv('reciprocal_h-m.txt', sep=" ", header=None)
data=data[~data[0].isin(gene_names)]
data=data[~data[1].isin(gene_names)]
data=data[~data[0].isin(gene_names_mac)]
data=data[~data[1].isin(gene_names_mac)]

data.to_csv('~/conserved_gene_order_demo/macaca_reference/Reciprocal_h-m.txt', header=None, index=None, sep=' ')

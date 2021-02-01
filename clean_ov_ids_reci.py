import ast
import pandas as pd
import os

def clean_reciprocal(path1, path2, overlapped1, overlapped2, reciporcal, path_save):

    # Open the file with human proteins having overlapped ids- isoforms-.

    os.chdir(os.path.expanduser(path1))

    with open(overlapped1) as file:
        gene_names = ast.literal_eval(file.read())

    os.chdir(os.path.expanduser(path2))

    with open(overlapped2) as file:
        gene_names_mac = ast.literal_eval(file.read())

    # Open human-macaca blast results with pandas and delete those rows
    # that contain the overlapped proteins.

    data = pd.read_csv(reciporcal, sep=" ", header=None)
    data=data[~data[0].isin(gene_names)]
    data=data[~data[1].isin(gene_names)]
    data=data[~data[0].isin(gene_names_mac)]
    data=data[~data[1].isin(gene_names_mac)]

    data.to_csv(path_save, header=None, index=None, sep=' ')

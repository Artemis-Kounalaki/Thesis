import ast
import pandas as pd
import os

def clean(path,ids,results,path_save):

        # Open the file with proteins having overlapped ids.

        os.chdir(os.path.expanduser(path))
        with open(ids) as file:
            gene_names = ast.literal_eval(file.read())


        # Open human blast results with pandas and delete those rows
        # that contain the overlapped proteins.

        data = pd.read_csv(results, sep="\t", header=None)
        data=data[~data[0].isin(gene_names)]
        data=data[~data[1].isin(gene_names)]
        data.to_csv(path_save, header=None, index=None, sep=' ', mode='a')


        # HOW MANY UNIQUE IDS NOW REMAIN.

        #column_values = data[[0, 1]].values.ravel()
        #unique_values =  pd.unique(column_values)
        #print(len(unique_values))

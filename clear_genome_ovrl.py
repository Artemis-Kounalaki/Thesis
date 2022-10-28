import ast
import os

# After finfing the isoforms, with this function we clear genome from isoforms.

def clear_genome(path,reference,overlapped,save):
    os.chdir(os.path.expanduser(path))
    with open(reference) as f:
        list2 = f.readlines()

    with open(overlapped) as file:
        gene_names = ast.literal_eval(file.read())

    for line in range(len(list2)):
        if any(gene in list2[line] for gene in gene_names):
            list2[line]=' '
            if line == len(list2)-1:
                break
            list2[line+1]=' '

    with open(save, 'a') as outfile:
        for sublist in list2:
           outfile.write('{}\n'.format(sublist))

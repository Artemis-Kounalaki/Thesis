import os
import numpy as np
from collections import defaultdict
from itertools import chain

#os.chdir(os.path.expanduser('~/conserved_gene_order/macaca_reference'))

# Make dcitionary with keys the human proteins and values the 'similar' macaca proteins.

data = np.loadtxt('reciprocal_h-m.txt', dtype=str)
data2=np.sort(data[:,0:2],axis=1)

Diction = defaultdict(list)
for element in data2:
    Diction[element[1]].append(element[0])
new_Dic = {a:list(set(b)) for a, b in Diction.items()}


# Make dictionary with human proteins and their paralogs.
#os.chdir(os.path.expanduser('~/conserved_gene_order/human_reference'))
datahum=np.loadtxt('clean_blast_results.txt', dtype=str)

Dictionhum1 = defaultdict(list)
for row in datahum:
    Dictionhum1[row[0]].append(row[1])

Dictionhum2 = defaultdict(list)
for row in datahum:
    Dictionhum2[row[1]].append(row[0])

Dictionhum = defaultdict(list)
for k, v in chain(Dictionhum1.items(), Dictionhum2.items()):
    Dictionhum[k].append(', '.join(v))

Diction_hum = {a:list(set(b)) for a, b in Dictionhum.items()}


# Combine the two dictionaries : keys = human proteins values = 1) macaca proteins, 2) human paralogs

d_comb = {key:[new_Dic[key], Diction_hum[key]] for key in new_Dic}

#Open overlapped ids in pairs.

from collections import defaultdict
import numpy as np
import os
import ast
import itertools

os.chdir(os.path.expanduser("~/conserved_gene_order/human_reference"))

with open("overlapped_ids.txt") as file:
    mylist = ast.literal_eval(file.read())
    gene_names=list(itertools.chain.from_iterable(mylist))
one=gene_names[::2]
two=gene_names[1::2]


#Make a new list with pairs sorted.

new=[list(t) for t in zip(one, two)]
new1=sorted(new)
new1=list(itertools.chain.from_iterable(new))
new=[x+y for x,y in zip(new1[0::2], new1[1::2])]


#Make a dictionary with keys as pairs.

Dict1=defaultdict(list)
for row in new:
    Dict1[row].append(0)


#Read blast results - make numpy array sorted

data = np.loadtxt('results_human.txt', dtype=str)
data2=np.sort(data[:,0:2],axis=1)


#Merge the two transcript columns.

c=np.core.defchararray.add(data2[:,0],data2[:,1])
data=np.column_stack((c,data))


#Make a dictionary with keys be the sorted merged transcript names.

dict2=defaultdict(list)
for i in data:
    dict2[i[0]].append(list(i[1:]))


#Keep a dictionary with keys only appear
#in dict2 and not in Dict1 (clean overlapped transcripts)

keep={k:v for k,v in dict2.items() if k not in Dict1}
unkeep=list(itertools.chain.from_iterable(keep.values()))


#Make a txt file with the transformed blast results.
with open("clean_blast_res.txt", 'w') as out:
    for item in unkeep:
        out.write(str(item)+'\n')
out.close()

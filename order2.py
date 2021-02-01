import os
import numpy as np
from collections import defaultdict
import re
import itertools
import pandas as pd
import ast

os.chdir(os.path.expanduser('~/conserved_gene_order/mus_reference'))



# Make dcitionary with keys the human proteins and values the 'similar' macaca proteins.

data = np.loadtxt('Reciprocal_h-mus.txt', dtype=str)
data2=np.sort(data[:,0:2],axis=1)

Diction = defaultdict(list)
for element in data2:
    Diction[element[1]].append(element[0])

new_Dic = {a:list(set(b)) for a, b in Diction.items()}

pr_hum=[]
pr_mac=[]
for key, value in new_Dic.items():
    pr_hum.append(key)
    pr_mac.append(value)
print(len(pr_mac))
print(pr_hum[0:2])
print(pr_mac[0:2])


os.chdir(os.path.expanduser('~/conserved_gene_order/human_reference'))
print('START HUM')
org='GRCh38'
df_human = pd.DataFrame(columns=['ID','Chrom','Start','End'])
counter=-1
with open('ids.txt') as f:
    for i in f:
        counter+=1
        name = re.findall('>+(\w+.\d+)',i)[0]
        chromosome = re.findall('%s:+(\w+):'% org,i)[0]
        start = re.findall('%s:+\w+:(\d+):' % org,i)[0]
        end = re.findall('%s:+\w+:\d+:(\d+):'% org,i)[0]
        df_human.loc[counter, 'ID'] = name
        df_human.loc[counter, 'Chrom'] = chromosome
        df_human.loc[counter, 'Start'] = start
        df_human.loc[counter, 'End'] = end

# Sort by chromosome and start

df_human.Chrom = df_human.Chrom.astype(str)
df_human.Start = df_human.Start.astype(int)
df_human = df_human.sort_values(by=['Chrom','Start'])
df_human = df_human.reset_index(drop=True)
#df_human=df_human[~df_human.Chrom.str.contains("MT")]
#df_human = df_human.reset_index(drop=True)


with open("overlapped_ids.txt") as file:
    gene_names = ast.literal_eval(file.read())

df_human=df_human[~df_human.ID.isin(gene_names)]
df_human = df_human.reset_index(drop=True)


os.chdir(os.path.expanduser('~/conserved_gene_order/mus_reference'))
print('START MUS')
org='GRCm38'
df_macaca = pd.DataFrame(columns=['ID','Chrom','Start','End'])
counter=-1
with open('mus_ids.txt') as f:
    for i in f:
        counter+=1
        name = re.findall('>+(\w+.\d+)',i)[0]
        chromosome = re.findall('%s:+(\w+):'% org,i)[0]
        start = re.findall('%s:+\w+:(\d+):' % org,i)[0]
        end = re.findall('%s:+\w+:\d+:(\d+):'% org,i)[0]
        df_macaca.loc[counter, 'ID'] = name
        df_macaca.loc[counter, 'Chrom'] = chromosome
        df_macaca.loc[counter, 'Start'] = start
        df_macaca.loc[counter, 'End'] = end

# Sort by chromosome and start

df_macaca.Chrom = df_macaca.Chrom.astype(str)
df_macaca.Start = df_macaca.Start.astype(int)
df_macaca = df_macaca.sort_values(by=['Chrom','Start'])
df_macaca = df_macaca.reset_index(drop=True)

with open("overlapped_ids_mus.txt") as file:
    gene_names = ast.literal_eval(file.read())

df_macaca=df_macaca[~df_macaca.ID.isin(gene_names)]
df_macaca = df_macaca.reset_index(drop=True)
print('OK')


# It's time to find which pairs are CGO and nCGO

cCGO=-1
cnCGO=-1
df_CGO= pd.DataFrame(columns=['Prot_human','Prot_mus'])
df_nCGO= pd.DataFrame(columns=['Prot_human','Prot_mus'])
seen=[]
counter=range(len(pr_mac)-1)
for i in counter:
    if i%1000 == 0:
        print('1000')
    for p in pr_mac[i]:
        pr_h=pr_hum[i]
        pr_m=p
        if [pr_h,pr_m] in seen:
            break
        seen.append([pr_h,pr_m])
        ind_h=df_human[df_human['ID']==pr_h].index.values[0]
        ind_m=df_macaca[df_macaca['ID']==pr_m].index.values[0]
        if ind_h==0 and ind_m==0:
            listA1 = [i for i, el in enumerate(pr_hum) if df_human.loc[ind_h+1,'ID'] in el]
            listB1 = [i for i, el in enumerate(pr_mac) if df_macaca.loc[ind_m+1,'ID'] in el]
            if any(x in listA1 for x in listB1):
                   cCGO+=1
                   df_CGO.loc[cCGO, 'Prot_human'] = pr_h
                   df_CGO.loc[cCGO, 'Prot_mus'] = pr_m
            else:
                cnCGO+=1
                df_nCGO.loc[cnCGO, 'Prot_human'] = pr_h
                df_nCGO.loc[cnCGO, 'Prot_mus'] = pr_m
        elif ind_h==0 and ind_m!=0 and ind_m!=len(df_macaca)-1:
            listA1 = [i for i, el in enumerate(pr_hum) if df_human.loc[ind_h+1,'ID'] in el]
            listB_1 = [i for i, el in enumerate(pr_mac) if df_macaca.loc[ind_m-1,'ID'] in el]
            listB1 = [i for i, el in enumerate(pr_mac) if df_macaca.loc[ind_m+1,'ID'] in el]
            if (any(x in listA1 for x in listB_1) or any(x in listA1 for x in listB1)) :
                   cCGO+=1
                   df_CGO.loc[cCGO, 'Prot_human'] = pr_h
                   df_CGO.loc[cCGO, 'Prot_mus'] = pr_m
            else:
                cnCGO+=1
                df_nCGO.loc[cnCGO, 'Prot_human'] = pr_h
                df_nCGO.loc[cnCGO, 'Prot_mus'] = pr_m
        elif ind_h==0 and ind_m==len(df_macaca)-1:
            listA1 = [i for i, el in enumerate(pr_hum) if df_human.loc[ind_h+1,'ID'] in el]
            listB_1 = [i for i, el in enumerate(pr_mac) if df_macaca.loc[ind_m-1,'ID'] in el]
            if any(x in listA1 for x in listB_1):
                cCGO+=1
                df_CGO.loc[cCGO, 'Prot_human'] = pr_h
                df_CGO.loc[cCGO, 'Prot_mus'] = pr_m
            else:
                cnCGO+=1
                df_nCGO.loc[cnCGO, 'Prot_human'] = pr_h
                df_nCGO.loc[cnCGO, 'Prot_mus'] = pr_m
        elif ind_h!=0 and ind_h!=len(df_human)-1 and ind_m==0:
            listA_1 =[i for i, el in enumerate(pr_hum) if df_human.loc[ind_h-1,'ID'] in el]
            listA1 = [i for i, el in enumerate(pr_hum) if df_human.loc[ind_h+1,'ID'] in el]
            listB1 = [i for i, el in enumerate(pr_mac) if df_macaca.loc[ind_m+1,'ID'] in el]
            if (any(x in listA1 for x in listB1) or any(x in listA_1 for x in listB1)):
                   cCGO+=1
                   df_CGO.loc[cCGO, 'Prot_human'] = pr_h
                   df_CGO.loc[cCGO, 'Prot_mus'] = pr_m
            else:
                cnCGO+=1
                df_nCGO.loc[cnCGO, 'Prot_human'] = pr_h
                df_nCGO.loc[cnCGO, 'Prot_mus'] = pr_m
        elif ind_h==len(df_human)-1 and ind_m==0:
            listA_1 = [i for i, el in enumerate(pr_hum) if df_human.loc[ind_h-1,'ID'] in el]
            listB1 = [i for i, el in enumerate(pr_mac) if df_macaca.loc[ind_m+1,'ID'] in el]
            if any(x in listA_1 for x in listB1) :
                   cCGO+=1
                   df_CGO.loc[cCGO, 'Prot_human'] = pr_h
                   df_CGO.loc[cCGO, 'Prot_mus'] = pr_m
            else:
                cnCGO+=1
                df_nCGO.loc[cnCGO, 'Prot_human'] = pr_h
                df_nCGO.loc[cnCGO, 'Prot_mus'] = pr_m
        elif ind_h==len(df_human)-1 and ind_m==len(df_macaca)-1:
            listA_1 = [i for i, el in enumerate(pr_hum) if df_human.loc[ind_h-1,'ID'] in el]
            listB_1 = [i for i, el in enumerate(pr_mac) if df_macaca.loc[ind_m-1,'ID'] in el]
            if any(x in listA_1 for x in listB_1):
                cCGO+=1
                df_CGO.loc[cCGO, 'Prot_human'] = pr_h
                df_CGO.loc[cCGO, 'Prot_mus'] = pr_m
            else:
                cnCGO+=1
                df_nCGO.loc[cnCGO, 'Prot_human'] = pr_h
                df_nCGO.loc[cnCGO, 'Prot_mus'] = pr_m
        elif ind_h!=0 and ind_h!=len(df_human)-1 and ind_m!=0 and ind_m!=len(df_macaca)-1:
            listA_1 =[i for i, el in enumerate(pr_hum) if df_human.loc[ind_h-1,'ID'] in el]
            listA1 = [i for i, el in enumerate(pr_hum) if df_human.loc[ind_h+1,'ID'] in el]
            listB_1 = [i for i, el in enumerate(pr_mac) if df_macaca.loc[ind_m-1,'ID'] in el]
            listB1 = [i for i, el in enumerate(pr_mac) if df_macaca.loc[ind_m+1,'ID'] in el]

            if (any(x in listA_1 for x in listB_1) or any(x in listA_1 for x in listB1) or any(x in listA1 for x in listB_1) or any(x in listA1 for x in listB1)):
                cCGO+=1
                df_CGO.loc[cCGO, 'Prot_human'] = pr_h
                df_CGO.loc[cCGO, 'Prot_mus'] = pr_m
            else:
                cnCGO+=1
                df_nCGO.loc[cnCGO, 'Prot_human'] = pr_h
                df_nCGO.loc[cnCGO, 'Prot_mus'] = pr_m

df_CGO.to_csv('~/conserved_gene_order/mus_reference/CGO.txt', header=None, index=None, sep=' ', mode='a')
df_nCGO.to_csv('~/conserved_gene_order/mus_reference/nCGO.txt', header=None, index=None, sep=' ', mode='a')

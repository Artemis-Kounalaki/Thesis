import os
import numpy as np
from collections import defaultdict
import re
import itertools
import pandas as pd

os.chdir(os.path.expanduser('~/conserved_gene_order/macaca_reference'))


# Make dcitionary with keys the human proteins and values the 'similar' macaca proteins.

data = np.loadtxt('Reciprocal_h-m.txt', dtype=str)
data2=np.sort(data[:,0:2],axis=1)

Diction = defaultdict(list)
for element in data2:
    Diction[element[1]].append(element[0])

new_Dic = {a:list(set(b)) for a, b in Diction.items()}



# Make dictionary with human proteins and their paralogs.
os.chdir(os.path.expanduser('~/conserved_gene_order/human_reference'))


datahum=np.loadtxt('clean_blast_results_human.txt', dtype=str)
reciprocal_hum=list(data2[:,1])
df=pd.DataFrame(data=datahum[0:,0:])
df=df[df[[0, 1]].isin(reciprocal_hum).all(1)]
datahum=df.to_numpy()

Dictionhum1 = defaultdict(list)
for row in datahum:
    Dictionhum1[row[0]].append(row[1])

Dictionhum2 = defaultdict(list)
for row in datahum:
    Dictionhum2[row[1]].append(row[0])

Dictionhum = Dictionhum1.copy()
Dictionhum.update(Dictionhum2)


Diction_hum = {a:list(set(b)) for a, b in Dictionhum.items()}


# Combine the two dictionaries : keys = human proteins values = 1) macaca proteins, 2) human paralogs

d_comb = {key:[new_Dic[key], Diction_hum[key]] for key in new_Dic}

s=np.hstack((datahum[:,0],datahum[:,1]))


# Make one list with paralog human proteins and one other with macaca proteins with similarity.

pr_hum=[]
pr_mac=[]

for key, value in d_comb.items():
    pr_hum.append(list(set([key]+value[1])))
    pr_mac.append(value[0])


df_human = pd.DataFrame(columns=['ID','Chrom','Start','End'])
f = open('ids.txt', 'r')
content = f.read()
f.close()
counter=-1
for i in list(set(list(itertools.chain.from_iterable(pr_hum)))):
    counter+=1
    line = re.findall(r"^>%s.*GRCh38:+(\w+)+:(\d+):(\d+):" %i, content, re.MULTILINE)[0]
    chrom=line[0]
    start=line[1]
    end=line[2]
    df_human.loc[counter, 'ID'] = i
    df_human.loc[counter, 'Chrom'] = chrom
    df_human.loc[counter, 'Start'] = start
    df_human.loc[counter, 'End'] = end



df_macaca = pd.DataFrame(columns=['ID','Chrom','Start','End'])
f = open('macaca_ids.txt', 'r')
content = f.read()
f.close()
counter=-1
for i in list(set(list(itertools.chain.from_iterable(pr_mac)))):
    counter+=1
    line = re.findall(r"^>%s.*Mmul_10:+(\w+)+:(\d+):(\d+):" %i, content, re.MULTILINE)[0]
    chrom=line[0]
    start=line[1]
    end=line[2]
    df_macaca.loc[counter, 'ID'] = i
    df_macaca.loc[counter, 'Chrom'] = chrom
    df_macaca.loc[counter, 'Start'] = start
    df_macaca.loc[counter, 'End'] = end
print(df_macaca)

# Sort by chromosome and start

df_human.Chrom = df_human.Chrom.astype(str)
df_human.Start = df_human.Start.astype(int)
df_human = df_human.sort_values(by=['Chrom','Start'])
df_human = df_human.reset_index(drop=True)
df_macaca.Chrom = df_macaca.Chrom.astype(str)
df_macaca.Start = df_macaca.Start.astype(int)
df_macaca = df_macaca.sort_values(by=['Chrom','Start'])
df_macaca = df_macaca.reset_index(drop=True)


# It's time to find which paris are CGO and nCGO

cCGO=-1
cnCGO=-1
df_CGO= pd.DataFrame(columns=['Prot_human','Prot_macaca'])
df_nCGO= pd.DataFrame(columns=['Prot_human','Prot_macaca'])
seen=[]
for i in pr_hum:
    for x in pr_mac:
        for pr in i:
            for p in x:
                pr_h=pr
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
                           df_CGO.loc[cCGO, 'Prot_macaca'] = pr_m
                    else:
                        cnCGO+=1
                        df_nCGO.loc[cnCGO, 'Prot_human'] = pr_h
                        df_nCGO.loc[cnCGO, 'Prot_macaca'] = pr_m
                elif ind_h==0 and ind_m!=0 and ind_m!=len(df_macaca)-1:
                    listA1 = [i for i, el in enumerate(pr_hum) if df_human.loc[ind_h+1,'ID'] in el]
                    listB_1 = [i for i, el in enumerate(pr_mac) if df_macaca.loc[ind_m-1,'ID'] in el]
                    listB1 = [i for i, el in enumerate(pr_mac) if df_macaca.loc[ind_m+1,'ID'] in el]
                    if (any(x in listA1 for x in listB_1) or any(x in listA1 for x in listB1)) :
                           cCGO+=1
                           df_CGO.loc[cCGO, 'Prot_human'] = pr_h
                           df_CGO.loc[cCGO, 'Prot_macaca'] = pr_m
                    else:
                        cnCGO+=1
                        df_nCGO.loc[cnCGO, 'Prot_human'] = pr_h
                        df_nCGO.loc[cnCGO, 'Prot_macaca'] = pr_m
                elif ind_h==0 and ind_m==len(df_macaca)-1:
                    listA1 = [i for i, el in enumerate(pr_hum) if df_human.loc[ind_h+1,'ID'] in el]
                    listB_1 = [i for i, el in enumerate(pr_mac) if df_macaca.loc[ind_m-1,'ID'] in el]
                    if any(x in listA1 for x in listB_1):
                        cCGO+=1
                        df_CGO.loc[cCGO, 'Prot_human'] = pr_h
                        df_CGO.loc[cCGO, 'Prot_macaca'] = pr_m
                    else:
                        cnCGO+=1
                        df_nCGO.loc[cnCGO, 'Prot_human'] = pr_h
                        df_nCGO.loc[cnCGO, 'Prot_macaca'] = pr_m
                elif ind_h!=0 and ind_h!=len(df_human)-1 and ind_m==0:
                    listA_1 =[i for i, el in enumerate(pr_hum) if df_human.loc[ind_h-1,'ID'] in el]
                    listA1 = [i for i, el in enumerate(pr_hum) if df_human.loc[ind_h+1,'ID'] in el]
                    listB1 = [i for i, el in enumerate(pr_mac) if df_macaca.loc[ind_m+1,'ID'] in el]
                    if (any(x in listA1 for x in listB1) or any(x in listA_1 for x in listB1)):
                           cCGO+=1
                           df_CGO.loc[cCGO, 'Prot_human'] = pr_h
                           df_CGO.loc[cCGO, 'Prot_macaca'] = pr_m
                    else:
                        cnCGO+=1
                        df_nCGO.loc[cnCGO, 'Prot_human'] = pr_h
                        df_nCGO.loc[cnCGO, 'Prot_macaca'] = pr_m
                elif ind_h==len(df_human)-1 and ind_m==0:
                    listA_1 = [i for i, el in enumerate(pr_hum) if df_human.loc[ind_h-1,'ID'] in el]
                    listB1 = [i for i, el in enumerate(pr_mac) if df_macaca.loc[ind_m+1,'ID'] in el]
                    if any(x in listA_1 for x in listB1) :
                           cCGO+=1
                           df_CGO.loc[cCGO, 'Prot_human'] = pr_h
                           df_CGO.loc[cCGO, 'Prot_macaca'] = pr_m
                    else:
                        cnCGO+=1
                        df_nCGO.loc[cnCGO, 'Prot_human'] = pr_h
                        df_nCGO.loc[cnCGO, 'Prot_macaca'] = pr_m
                elif ind_h==len(df_human)-1 and ind_m==len(df_macaca)-1:
                    listA_1 = [i for i, el in enumerate(pr_hum) if df_human.loc[ind_h-1,'ID'] in el]
                    listB_1 = [i for i, el in enumerate(pr_mac) if df_macaca.loc[ind_m-1,'ID'] in el]
                    if any(x in listA_1 for x in listB_1):
                        cCGO+=1
                        df_CGO.loc[cCGO, 'Prot_human'] = pr_h
                        df_CGO.loc[cCGO, 'Prot_macaca'] = pr_m
                    else:
                        cnCGO+=1
                        df_nCGO.loc[cnCGO, 'Prot_human'] = pr_h
                        df_nCGO.loc[cnCGO, 'Prot_macaca'] = pr_m
                elif ind_h!=0 and ind_h!=len(df_human)-1 and ind_m!=0 and ind_m!=len(df_macaca)-1:
                    listA_1 =[i for i, el in enumerate(pr_hum) if df_human.loc[ind_h-1,'ID'] in el]
                    listA1 = [i for i, el in enumerate(pr_hum) if df_human.loc[ind_h+1,'ID'] in el]
                    listB_1 = [i for i, el in enumerate(pr_mac) if df_macaca.loc[ind_m-1,'ID'] in el]
                    listB1 = [i for i, el in enumerate(pr_mac) if df_macaca.loc[ind_m+1,'ID'] in el]

                    if (any(x in listA_1 for x in listB_1) or any(x in listA_1 for x in listB1) or any(x in listA1 for x in listB_1) or any(x in listA1 for x in listB1)):
                        cCGO+=1
                        df_CGO.loc[cCGO, 'Prot_human'] = pr_h
                        df_CGO.loc[cCGO, 'Prot_macaca'] = pr_m
                    else:
                        cnCGO+=1
                        df_nCGO.loc[cnCGO, 'Prot_human'] = pr_h
                        df_nCGO.loc[cnCGO, 'Prot_macaca'] = pr_m

import os
import numpy as np
from collections import defaultdict
import re
import pandas as pd
import ast

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

Dictionhum1 = defaultdict(list)
for row in datahum:
    Dictionhum1[row[0]].append(row[1])

Dictionhum2 = defaultdict(list)
for row in datahum:
    Dictionhum2[row[1]].append(row[0])

#Dictionhum = defaultdict(list)
#for k, v in chain(Dictionhum1.items(), Dictionhum2.items()):
    #Dictionhum[k].append(" '".join(v))

Dictionhum = Dictionhum1.copy()
Dictionhum.update(Dictionhum2)


Diction_hum = {a:list(set(b)) for a, b in Dictionhum.items()}


# Combine the two dictionaries : keys = human proteins values = 1) macaca proteins, 2) human paralogs

d_comb = {key:[new_Dic[key], Diction_hum[key]] for key in new_Dic}


# Make one list with paralog human proteins and one other with macaca proteins with similarity.
pr_hum=[]
pr_mac=[]

for key, value in d_comb.items():
    pr_hum.append(list(set([key]+value[1])))
    pr_mac.append(value[0])

print(len(pr_mac))
print('START HUM')

org='GRCh38'
df_human = pd.DataFrame(columns=['ID','Chrom','Start','End'])
#os.chdir(os.path.expanduser(path))
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
df_human=df_human[~df_human.Chrom.str.contains("MT")]
df_human = df_human.reset_index(drop=True)


with open("overlapped_ids.txt") as file:
    gene_names = ast.literal_eval(file.read())

df_human=df_human[~df_human.ID.isin(gene_names)]
df_human = df_human.reset_index(drop=True)



print('START MAC')
org='Mmul_10'
df_macaca = pd.DataFrame(columns=['ID','Chrom','Start','End'])
os.chdir(os.path.expanduser('~/conserved_gene_order/macaca_reference'))
counter=-1
with open('macaca_ids.txt') as f:
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
df_macaca=df_macaca[~df_macaca.Chrom.str.contains("MT")]
df_macaca = df_macaca.reset_index(drop=True)


with open("overlapped_ids_macaca.txt") as file:
    gene_names = ast.literal_eval(file.read())

df_macaca=df_macaca[~df_macaca.ID.isin(gene_names)]
df_macaca = df_macaca.reset_index(drop=True)
print('OK')

# Now I want to fix a function found when have 2 proteins the proteins included when CGO is present.

up_down=1
listhum=[]
listmac=[]
a1='ENSP00000374864.2'
b1='ENSMMUP00000079563.1'
listhum.append(a1)
listmac.append(b1)
index_h=df_human[df_human['ID']==a1].index.values[0]
index_m=df_macaca[df_macaca['ID']==b1].index.values[0]
print(index_m,index_h)
while True:
    if index_h!=0 and index_h!=len(df_human)-1 and index_m!=0 and index_m!=len(df_macaca)-1:
        listA_1 =[i for i, el in enumerate(pr_hum) if df_human.loc[index_h-up_down,'ID'] in el]
        listA1 = [i for i, el in enumerate(pr_hum) if df_human.loc[index_h+up_down,'ID'] in el]
        listB_1 = [i for i, el in enumerate(pr_mac) if df_macaca.loc[index_m-up_down,'ID'] in el]
        listB1 = [i for i, el in enumerate(pr_mac) if df_macaca.loc[index_m+up_down,'ID'] in el]
        if any(x in listA_1 for x in listB_1):
            listhum.insert(0,df_human.loc[index_h-up_down,'ID'])
            listmac.insert(0,df_macaca.loc[index_m-up_down,'ID'])
            up_down+=1
        elif any(x in listA_1 for x in listB1):
            listhum.insert(0,df_human.loc[index_h-up_down,'ID'])
            listmac.insert(len(listmac),df_macaca.loc[index_m+up_down,'ID'])
            up_down+=1
        elif any(x in listA1 for x in listB_1):
            listhum.insert(len(listhum),df_human.loc[index_h+up_down,'ID'])
            listmac.insert(0,df_macaca.loc[index_m-up_down,'ID'])
            up_down+=1
        elif any(x in listA1 for x in listB1):
            listhum.insert(len(listhum),df_human.loc[index_h+up_down,'ID'])
            listmac.insert(len(listmac),df_macaca.loc[index_m+up_down,'ID'])
            up_down+=1
        else:
            break
    elif index_h==0 and index_m==0:
        listA1 = [i for i, el in enumerate(pr_hum) if df_human.loc[index_h+up_down,'ID'] in el]
        listB1 = [i for i, el in enumerate(pr_mac) if df_macaca.loc[index_m+up_down,'ID'] in el]
        if any(x in listA1 for x in listB1):
            listhum.insert(len(listhum),df_human.loc[index_h+up_down,'ID'])
            listmac.insert(len(listmac),df_macaca.loc[index_m+up_down,'ID'])
            up_down+=1
        else:
            break
    elif index_h==0 and index_m!=0 and index_m!=len(df_macaca)-1:
        listA1 = [i for i, el in enumerate(pr_hum) if df_human.loc[index_h+up_down,'ID'] in el]
        listB_1 = [i for i, el in enumerate(pr_mac) if df_macaca.loc[index_m-up_down,'ID'] in el]
        listB1 = [i for i, el in enumerate(pr_mac) if df_macaca.loc[index_m+up_down,'ID'] in el]
        if any(x in listA1 for x in listB_1):
            listhum.insert(len(listhum),df_human.loc[index_h+up_down,'ID'])
            listmac.insert(0,df_macaca.loc[index_m-up_down,'ID'])
            up_down+=1
        elif any(x in listA1 for x in listB1) :
            listhum.insert(len(listhum),df_human.loc[index_h+up_down,'ID'])
            listmac.insert(len(listmac),df_macaca.loc[index_m+up_down,'ID'])
            up_down+=1
        else:
            break
    elif index_h==0 and index_m==len(df_macaca)-1:
        listA1 = [i for i, el in enumerate(pr_hum) if df_human.loc[index_h+up_down,'ID'] in el]
        listB_1 = [i for i, el in enumerate(pr_mac) if df_macaca.loc[index_m-up_down,'ID'] in el]
        if any(x in listA1 for x in listB_1):
            listhum.insert(len(listhum),df_human.loc[index_h+up_down,'ID'])
            listmac.insert(0,df_macaca.loc[index_m-up_down,'ID'])
            up_down+=1
        else:
            break
    elif index_h!=0 and index_h!=len(df_human)-1 and index_m==0:
        listA_1 =[i for i, el in enumerate(pr_hum) if df_human.loc[index_h-up_down,'ID'] in el]
        listA1 = [i for i, el in enumerate(pr_hum) if df_human.loc[index_h+up_down,'ID'] in el]
        listB1 = [i for i, el in enumerate(pr_mac) if df_macaca.loc[index_m+up_down,'ID'] in el]
        if any(x in listA1 for x in listB1):
            listhum.insert(len(listhum),df_human.loc[index_h+up_down,'ID'])
            listmac.insert(len(listmac),df_macaca.loc[index_m+up_down,'ID'])
            up_down+=1
        elif any(x in listA_1 for x in listB1):
            listhum.insert(0,df_human.loc[index_h-up_down,'ID'])
            listmac.insert(len(listmac),df_macaca.loc[index_m+up_down,'ID'])
            up_down+=1
        else:
            break
    elif index_h==len(df_human)-1 and index_m==0:
        listA_1 = [i for i, el in enumerate(pr_hum) if df_human.loc[index_h-up_down,'ID'] in el]
        listB1 = [i for i, el in enumerate(pr_mac) if df_macaca.loc[index_m+up_down,'ID'] in el]
        if any(x in listA_1 for x in listB1) :
            listhum.insert(0,df_human.loc[index_h-up_down,'ID'])
            listmac.insert(len(listmac),df_macaca.loc[index_m+up_down,'ID'])
            up_down+=1
        else:
            break
    elif index_h==len(df_human)-1 and index_m==len(df_macaca)-1:
        listA_1 = [i for i, el in enumerate(pr_hum) if df_human.loc[index_h-up_down,'ID'] in el]
        listB_1 = [i for i, el in enumerate(pr_mac) if df_macaca.loc[index_m-up_down,'ID'] in el]
        if any(x in listA_1 for x in listB_1):
            listhum.insert(0,df_human.loc[index_h-up_down,'ID'])
            listmac.insert(0,df_macaca.loc[index_m-up_down,'ID'])
            up_down+=1
        else:
            break

print(listmac)
print(listhum)
print(up_down)
print(len(listhum))
print(len(listmac))

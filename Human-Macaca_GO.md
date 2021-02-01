---
title: "Gene order conservation-sequel- (1)"
author: "Artemis Kounalaki"
---

**Now, its time to find ortholgs.** <br />
## 1st group: Homo sapiens - Macaca mulatta <br />
**1. Download macaca's protein whole genome sequences.**
<br />

```

mkdir $HOME/conserved_gene_order/macaca_reference
cd $HOME/conserved_gene_order/macaca_reference
wget 'ftp://ftp.ensembl.org/pub/release-100/fasta/macaca_mulatta/pep/Macaca_mulatta.Mmul_10.pep.all.fa.gz'
gunzip Macaca_mulatta.Mmul_10.pep.all.fa.gz


# Make sequences in one line

awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < Macaca_mulatta.Mmul_10.pep.all.fa > macaca_reference.fa

```


**2. Clean macaca's protein whole genome sequences.**

It's appropriate to clean file from sequences like scaffold regions, isoforms and unknown regions.<br/>

```

# scaffolds

awk '/scaffold/{n=2}; n {n--; next}; 1' < macaca_reference.fa > macaca_reference_cl.fa


# isoforms

awk '/isoform/{n=2}; n {n--; next}; 1' < macaca_reference_cl.fa > macaca_reference_cln.fa

# unkonown regions

awk '/QNVO/{n=2}; n {n--; next}; 1' < macaca_reference_cln.fa > macaca_reference_clean1.fa

# regions with L1 position

awk '/Mmul_10:ML1/{n=2}; n {n--; next}; 1' < macaca_reference_clean1.fa > macaca_reference_clean.fa



# remove useless files

rm macaca_reference.fa
rm macaca_reference_cl.fa
rm macaca_reference_cln.fa
rm macaca_reference_clean1.fa

```


**3. Extract protein ids and their location.**
<br/>

```
#Find ids
grep '>' macaca_reference_clean.fa> macaca_ids.txt

```


**4. Make Blast database and run blast.**

Run blast reciprocal to find macaca's proteins that are close enough to human proteins.
<br/>
```
# Make the blast database.

cd $HOME/conserved_gene_order/blast_db
makeblastdb -in $HOME/conserved_gene_order/macaca_reference/macaca_reference_clean.fa -dbtype prot -parse_seqids -out database_macaca
cd $HOME/conserved_gene_order/macaca_reference


# Run blastp. Human - Macaca

blastp -num_threads 16 -db $HOME/conserved_gene_order/blast_db/database_macaca -evalue 1e-10 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen" -qcov_hsp_perc 80 -query $HOME/conserved_gene_order/human_reference/human_reference_clean.fa >'results_human-macaca.txt'


# Run blastp. Macaca - Human

blastp -num_threads 16 -db $HOME/conserved_gene_order/blast_db/database_human -evalue 1e-10 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen" -qcov_hsp_perc 80 -query $HOME/conserved_gene_order/macaca_reference/macaca_reference_clean.fa >'results_macaca-human.txt'


# Run blastp. Macaca - Macaca

blastp -num_threads 16 -db $HOME/conserved_gene_order/blast_db/database_macaca -evalue 1e-10 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen" -qcov_hsp_perc 80 -query $HOME/conserved_gene_order/macaca_reference/macaca_reference_clean.fa >'results_macaca.txt'
cd $HOME
```


**5. Clean macaca blast results from isoforms.**
Clean blast results from isoforms with the same code we used to clean human blast results in Human_paralogs.
 <br/>
 This implementation is python3 script. <br/>
 1. Find the isoforms. _overlapped_ids.py_  <br/>
 2. Clean the results from the isoforms.  _clean_ov_ids.py_ <br/>

**6. Find reciprocal proteins(human-macaca & macaca-human).**
It's appropriate to exclude proteins that not match reciprocally. _reciprocal_mac-hum.py_ from _reciprocal.py_

```
from collections import defaultdict
import numpy as np
import os
import itertools

def reciprocal(path1, blast_results1,blast_results2, save1, save2) :

    #Read blast results - make numpy array sorted

    os.chdir(os.path.expanduser(path1))

    data1 = np.loadtxt(blast_results1, dtype=str)
    data2=np.sort(data1[:,0:2],axis=1)
    data11 =np.loadtxt(blast_results2, dtype=str)
    data22=np.sort(data11[:,0:2],axis=1)

    #Merge the two protein columns.

    c=np.core.defchararray.add(data2[:,0],data2[:,1])
    data1=np.column_stack((c,data1))

    c1=np.core.defchararray.add(data22[:,0],data22[:,1])
    data11=np.column_stack((c1,data11))

    #Make a dictionary with keys be the sorted merged transcript names.

    dict2=defaultdict(list)
    for i in data1:
        dict2[i[0]].append(list(i[1:]))

    Dict1=defaultdict(list)
    for i in data11:
        Dict1[i[0]].append(list(i[1:]))


    #Keep a dictionary with keys both appear
    #in dict2 and Dict1

    keep1={k:v for k,v in dict2.items() if k in Dict1}
    unkeep1=list(itertools.chain.from_iterable(keep1.values()))
    keep2={k:v for k,v in Dict1.items() if k in dict2}
    unkeep2=list(itertools.chain.from_iterable(keep2.values()))


    #Make a txt file with the transformed blast results.
    with open(save1, 'w') as out:
        for item in unkeep1:
            out.write(str(item)+'\n')
    out.close()
    with open(save2, 'w') as out:
        for item in unkeep2:
            out.write(str(item)+'\n')
    out.close()
 ```
 <br/>

Combine these two files. <br/>

```
cd $HOME/conserved_gene_order/macaca_reference
cat reciprocal_hum-mac*.txt >> reciprocal_hum-mac_all.txt
tr -d "['],"< reciprocal_hum-mac_all.txt > reciprocal_h-m.txt
rm reciprocal_hum-mac*.txt
cd $HOME
```

**7. Clean reciprocal results from human and macaca isoforms.**
Clean from overlapped ids found in macaca and in human. <br/>
_clean_ov_ids_h-mac.py_ from _clean_ov_ids_reci.py_ <br/>
```
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
```

**8. Find couples with conserved gene order.**<br/>
It's time to find which couples of proteins (human_protein - macaca_protein) have conserved gene order and which are not and save them at two different files.<br/>
Implementation code: **Python3** _order1.py_ <br/>

```
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


os.chdir(os.path.expanduser('~/conserved_gene_order/macaca_reference'))
print('START MAC')
org='Mmul_10'
df_macaca = pd.DataFrame(columns=['ID','Chrom','Start','End'])
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

with open("overlapped_ids_macaca.txt") as file:
    gene_names = ast.literal_eval(file.read())

df_macaca=df_macaca[~df_macaca.ID.isin(gene_names)]
df_macaca = df_macaca.reset_index(drop=True)
print('OK')


# It's time to find which pairs are CGO and nCGO

cCGO=-1
cnCGO=-1
df_CGO= pd.DataFrame(columns=['Prot_human','Prot_macaca'])
df_nCGO= pd.DataFrame(columns=['Prot_human','Prot_macaca'])
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

df_CGO.to_csv('~/conserved_gene_order/macaca_reference/CGO.txt', header=None, index=None, sep=' ', mode='a')
df_nCGO.to_csv('~/conserved_gene_order/macaca_reference/nCGO.txt', header=None, index=None, sep=' ', mode='a')

```
**9. Create a table with CGOs and nCGOs and their average percentage idenity in sequence level.**<br/>
Implementation code: **Python3** _identity_per.py_ <br/>
```
import pandas as pd
import numpy as np

os.chdir(os.path.expanduser('~/conserved_gene_order/macaca_reference'))
cgo = pd.read_csv('CGO.txt', sep=" ", header=None)
cgo.columns = ["Prot_human", "Prot_macaca"]
cgo = cgo.sort_values(by=['Prot_human'])
cgo = cgo.reset_index(drop=True)

ncgo = pd.read_csv('nCGO.txt', sep=" ", header=None)
ncgo.columns = ["Prot_human", "Prot_macaca"]
ncgo = ncgo.sort_values(by=['Prot_human'])
ncgo = ncgo.reset_index(drop=True)

reci = pd.read_csv('Reciprocal_h-m.txt', sep=" ", header=None)
reci.columns = ["Prot_human", "Prot_macaca","P.identity", 'Length', '4', 'Gap', '6', '7', '8', '9', '10', '11', 'Hum_len', 'Mac_len']


Cgo = pd.merge(reci, cgo, on=['Prot_human','Prot_macaca'], how='left', indicator='Exist')
Cgo['Exist'] = np.where(Cgo.Exist == 'both', True, False)

Cgo=Cgo.loc[Cgo['Exist'] == True]
Cgo=Cgo.drop(columns=['4', '6', '7', '8', '9', '10', '11'])
Cgo = Cgo.drop_duplicates(subset = ["Prot_human",'Prot_macaca'])
Cgo = Cgo.reset_index(drop=True)
Cgo['Conserved']='1'
nCgo = pd.merge(reci, ncgo, on=['Prot_human','Prot_macaca'], how='left', indicator='Exist')
nCgo['Exist'] = np.where(nCgo.Exist == 'both', True, False)

nCgo=nCgo.loc[nCgo['Exist'] == True]
nCgo=nCgo.drop(columns=[ '4','6', '7', '8', '9', '10', '11'])
nCgo = nCgo.drop_duplicates(subset = ["Prot_human",'Prot_macaca'])
nCgo = nCgo.reset_index(drop=True)
nCgo['Conserved']='0'

df_al=pd.concat([Cgo, nCgo], ignore_index=True)
df_al=df_al.drop(columns=['Exist'])
df_al = df_al.reset_index(drop=True)


print(Cgo['P.identity'].sum()/len(Cgo))
print(nCgo['P.identity'].sum()/len(nCgo))

df_al.to_csv('al.txt', sep='\t')
```

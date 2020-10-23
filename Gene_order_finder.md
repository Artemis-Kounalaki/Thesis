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


#make sequences in one line

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

awk '/QNVO/{n=2}; n {n--; next}; 1' < macaca_reference_cln.fa > macaca_reference_clean.fa


# remove useless files

rm macaca_reference.fa
rm macaca_reference_cl.fa
rm macaca_reference_cln.fa

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
cd $HOME
```


**5. Clear the results.**
In order to find human proteins that are very close to macaca proteins, it's appropriate to keep only those proteins that are found reciprocally.
 <br/>
 This implementation is python3 script. <br/>

 ```
 from collections import defaultdict
 import numpy as np
 import os
 import itertools



 #Read blast results - make numpy array sorted

 os.chdir(os.path.expanduser('~/conserved_gene_order/macaca_reference'))

 data1 = np.loadtxt('results_human-macaca.txt', dtype=str)
 data2=np.sort(data1[:,0:2],axis=1)
 data11 =np.loadtxt('results_macaca-human.txt', dtype=str)
 data22=np.sort(data11[:,0:2],axis=1)

 #Merge the two transcript columns.

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
 with open("reciprocal_hum-mac1.txt", 'w') as out:
     for item in unkeep1:
         out.write(str(item)+'\n')
 out.close()
 with open("reciprocal_hum-mac2.txt", 'w') as out:
     for item in unkeep2:
         out.write(str(item)+'\n')
 out.close()
 ```
 <br/>

Combine these two files. <br/>
```
cat reciprocal_hum-mac*.txt >> reciprocal_hum-mac_all.txt
tr -d "['],"< reciprocal_hum-mac_all.txt > reciprocal_h-m.txt
rm reciprocal_hum-mac.*.txt
```

**6.Grouping proteins.**

Make groups: human proteins that are similar enough with macaca's proteins.
<br/>
Implementation code: Python3 <br/>

```
import numpy as np
from collections import defaultdict
from itertools import chain


# Make dcitionary with keys the human proteins and values the 'similar' macaca proteins.

data = np.loadtxt('reciprocal_h-m.txt', dtype=str)
data2=np.sort(data[:,0:2],axis=1)

Diction = defaultdict(list)
for element in data2:
    Diction[element[1]].append(element[0])
new_Dic = {a:list(set(b)) for a, b in Diction.items()}


# Make dictionary with human proteins and their homologs.

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


# Combine the two dictionaries : keys = human proteins values = 1) macaca proteins, 2) human homologs

d_comb = {key:[new_Dic[key], Diction_hum[key]] for key in new_Dic}
```

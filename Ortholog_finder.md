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

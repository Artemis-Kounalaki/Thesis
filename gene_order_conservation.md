---
title: "Gene order conservation"
author: "Artemis Kounalaki"
output: html_document
---

<br />
**1. Download the data**

In order to finding duplicated regions in human, protein sequences are downloaded from Ensembl db.
<br />
```
# Make a file for the whole process.

mkdir $HOME/conserved_gene_order
cd $HOME/conserved_gene_order


# Make a file inside with the human protein sequences from Ensembl.

mkdir human_reference
cd human_reference
wget 'ftp://ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz'
gunzip Homo_sapiens.GRCh38.pep.all.fa.gz

```

**2. This is the code for cleaning the fasta file of human protein sequences.**

Find scaffold regions - unknown regions-, isoforms and haplotypic regions.
Make a file with these lines.
Remove these lines from the inital file with human protein sequences.

<br />

```
#scaffolds
cd $HOME/conserved_gene_order/human_reference
grep -A 1 -w 'scaffold' Homo_sapiens.GRCh38.pep.all.fa> scaf.sed
sed '/--/d' ./scaf.sed> scaff.sed
rm scaf.sed
grep -Fvx -f scaff.sed Homo_sapiens.GRCh38.pep.all.fa >reference.fa


# isoforms
grep -A 1 -w 'isoform' reference.fa> isoform.sed
grep -Fvx -f isoform.sed reference.fa > human_reference.fa


# haplotypic regions
awk '/CHR_/{n=2}; n {n--; next}; 1' < human_reference.fa > human_reference_cl.fa

#awk '/CHR_HSCHR/{n=2}; n {n--; next}; 1' < human_reference.fa > human_reference_cl.fa
#grep -A1 'CHR_HSCHR' human_reference.fa> chr.sed
#sed '/--/d' ./chr.sed> chrr.sed
#rm chr.sed
#grep -Fvx -f chrr.sed human_reference.fa >human_reference_cl.fa


# remove useless files
rm isoform.sed
rm reference.fa
rm human_reference.fa
rm scaff.sed

```

**3. Install Blast**

Install Blast and create a database with the reference sequences.
<br />

```
# Blast installation from ncbi

cd $HOME/conserved_gene_order
#wget https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi-blast-2.10.1+-x64-macosx.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi-blast-2.10.1+-x64-linux.tar.gz
tar -xvzf ncbi-blast-2.10.1+-x64-linux.tar.gz
#tar -xvzf ncbi-blast-2.10.1+-x64-macosx.tar.gz
rm ncbi-blast-2.10.1+-x64-linux.tar.gz
cd $HOME/conserved_gene_order/ncbi-blast-2.10.1+/bin
export PATH=$PATH:$pwd


# Database

cd ..
mkdir $HOME/conserved_gene_order/blast_db
export BLASTDB=$HOME/conserved_gene_order/blast_db
set BLASTDB=$HOME/conserved_gene_order/blast_db
cd $HOME/conserved_gene_order/blast_db


# Make the blast database.

makeblastdb -in $HOME/conserved_gene_order/human_reference/human_reference_cl.fa -dbtype prot -out database_human


# Run blastp.

blastp -num_threads 16 -db $HOME/conserved_gene_order/blast_db/database_human -evalue 1e-10 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen" -qcov_hsp_p 80 -query $HOME/conserved_gene_order/human_reference/human_reference.fa >'results_human.txt'

```

**4. Clean Blast results**

In order to avoid finding overlapped genes, it is appropriate to exclude those genes whose position in chromosome is inside another gene's position.
<br />
So, we get those ids that are used in Blastp.

```
cd $HOME/conserved_gene_order/reference_human
grep '>' human_reference_cl.fa> ids.txt

```
<br />

#Create a list with those genes that are overlapped in pairs. The lists include: gene1&2  name, chrosmosome, start, end.
#Python3 code is given.

```
import os
import re
myis=[]
with open("ids.txt") as f:
    for i in f:
        name = re.findall('>+(\w+.\d+)',i)[0]
        chromosome = re.findall('GRCh38:+(\w+):',i)[0]
        start = re.findall('GRCh38:+\w+:(\d+):',i)[0]
        end = re.findall('GRCh38:+\w+:\d+:(\d+):',i)[0]
        myis.append([name,chromosome ,start, end])


os.chdir("$HOME/conserved_gene_order/blast_db")
id1=[]
id2=[]
with open("results_human.txt") as file:
    for line in file:
        id1.append(line.split('\t')[0])
        id2.append(line.split('\t')[1])
pairs=list(zip(id1,id2))


njobs=-1
# in which list id is included
overlapped_pairs=[]
for pair in pairs:
    a= [i for i, el in enumerate(myis) if pair[0] in el][0]
    b= [i for i, el in enumerate(myis) if pair[1] in el][0]

    ch1=myis[a][1]
    st1=int(myis[a][2])
    end1=int(myis[a][3])
    range1=range(st1,end1)

    ch2=myis[b][1]
    st2=int(myis[b][2])
    end2=int(myis[b][3])
    range2=range(st2,end2)


    if ch1==ch2 and ((st1>st2 and st1<end2) or (st2>st1 and st2<end1)):
        overlapped_pairs.append([pair[0],pair[1]])

with open("overlapped_ids.txt", "w") as output:
    output.write(str(overlapped_pairs))
```

Now I have all trnascript ids with their location and the blast results.

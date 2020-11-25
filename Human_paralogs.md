---
title: "Gene order conservation"
author: "Artemis Kounalaki"
---

**1. Download the data**

In order to find duplicated regions in human, protein sequences are downloaded from Ensembl db.
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

Make multiline sequences into one line sequence each. <br />
Find scaffold regions - unknown regions-, isoforms and haplotypic regions. <br />
Remove these lines from the initial file with human protein sequences.
<br />

```
# Make sequences in one line

awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < Homo_sapiens.GRCh38.pep.all.fa > human_reference.fa


# scaffolds

awk '/scaffold/{n=2}; n {n--; next}; 1' < human_reference.fa >
human_reference_cl.fa


# isoforms

awk '/isoform/{n=2}; n {n--; next}; 1' < human_reference_cl.fa > human_reference_cln.fa


# haplotypic regions

awk '/CHR_/{n=2}; n {n--; next}; 1' < human_reference_cln.fa > human_reference_clean.fa


# remove useless files

rm human_reference.fa
rm human_reference_cl.fa
rm human_reference_cln.fa

```

**3. Install Blast**

Install Blast and create a database with the reference sequences.
<br />

```
# Blast installation from ncbi

cd $HOME/conserved_gene_order
#wget https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi-blast-2.11.0+-x64-macosx.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi-blast-2.11.0+-x64-linux.tar.gz
tar -xvzf ncbi-blast-2.11.0+-x64-linux.tar.gz
#tar -xvzf ncbi-blast-2.11.0+-x64-macosx.tar.gz
rm ncbi-blast-2.11.0+-x64-linux.tar.gz
cd $HOME/conserved_gene_order/ncbi-blast-2.11.0+/bin
export PATH=$PATH:$pwd


# Database

cd ..
mkdir $HOME/conserved_gene_order/blast_db
export BLASTDB=$HOME/conserved_gene_order/blast_db
set BLASTDB=$HOME/conserved_gene_order/blast_db


# Make the blast database.

cd $HOME/conserved_gene_order/blast_db
makeblastdb -in $HOME/conserved_gene_order/human_reference/human_reference_clean.fa -dbtype prot -parse_seqids  -out database_human


```
**3. Run Blast and find human paralogs.**

Run blastp.
<br />

```
cd $HOME/conserved_gene_order/human_reference
blastp -num_threads 16 -db $HOME/conserved_gene_order/blast_db/database_human -evalue 1e-10 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen"  -qcov_hsp_perc 80 -query $HOME/conserved_gene_order/human_reference/human_reference_clean.fa >'results_human.txt'

```

**4. Clean Blast results**

In order to avoid finding overlapped genes, it is appropriate to exclude those genes whose position in chromosome is inside another gene's position. <br />
So, we get those ids that are used in Blastp. <br />

```
cd $HOME/conserved_gene_order/human_reference
grep '>' human_reference_clean.fa> ids.txt
cd $HOME

```
Check every pair of blast results if proteins are overlapped. <br />
The lists include: gene1&2  name, chromosome, start, end. <br />
Save in a list the smallest protein of every overlapping pair in order to retain only the biggest protein and not the isoforms. <br />
**Python3** code is given. _Overlapped_ids.py_ <br />

```
import os
import re


def overlap(path,ids,org,blast_results,save_file):
    os.chdir(os.path.expanduser(path))
    myis=[]
    with open(ids) as f:
        for i in f:
            name = re.findall('>+(\w+.\d+)',i)[0]
            chromosome = re.findall('%s:+(\w+):'% org,i)[0]
            start = re.findall('%s:+\w+:(\d+):' % org,i)[0]
            end = re.findall('%s:+\w+:\d+:(\d+):'% org,i)[0]
            myis.append([name,chromosome ,start, end])


    id1=[]
    id2=[]
    with open(blast_results) as file:
        for line in file:
            id1.append(line.split('\t')[0])
            id2.append(line.split('\t')[1])
    pairs=list(zip(id1,id2))


    # in which list id is included

    overlapped_pr=[]
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
            if (end1-st1)>(end2-st2):
                overlapped_pr.append(pair[0])
            elif (end2-st2)>(end1-st1):
                overlapped_pr.append(pair[1])
            #overlapped_pairs.append([pair[0],pair[1]])

    with open(save_file, "w") as output:
        output.write(str(overlapped_pr))


```

Now I have all protein ids with their location and the blast results.<br />
It's time to clean the blast results from overlapped proteins (isoforms) in **Python3**. _Clean_ov_ids.py_ <br />

```
import ast
import pandas as pd
import os

def clean(path,ids,results,path_save):

        # Open the file with proteins having overlapped ids.

        os.chdir(os.path.expanduser(path))
        with open(ids) as file:
            gene_names = ast.literal_eval(file.read())


        # Open human blast results with pandas and delete those rows
        # that contain the overlapped proteins.

        data = pd.read_csv(results, sep="\t", header=None)
        data=data[~data[0].isin(gene_names)]
        data=data[~data[1].isin(gene_names)]
        data.to_csv(path_save, header=None, index=None, sep=' ', mode='a')


        # HOW MANY UNIQUE IDS NOW REMAIN.

        #column_values = data[[0, 1]].values.ravel()
        #unique_values =  pd.unique(column_values)
        #print(len(unique_values))

```

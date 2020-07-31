# Thesis
Gene order conservation
---
title: "Gene order conservation"
author: "Artemis Kounalaki"
output: html_document
---

<br />
**1. Download the data.**

In order to finding duplicated regions in human, protein sequences are downloaded from Ensembl db.
<br />
```{r eval=FALSE}
# Make a file for the whole process.

mkdir $HOME/conserved_gene_order
cd $HOME/conserved_gene_order


# Make a file inside with the human protein sequences from Ensembl.

mkdir human_reference
cd human_reference
wget 'ftp://ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz'
gunzip Homo_sapiens.GRCh38.pep.all.fa.gz

```
<br />
**2. This is the code for cleaning the fasta file of human protein sequences.**

Find scaffold regions - unknown regions-, isoforms and haplotypic regions.
Make a file with these lines.
Remove these lines from the inital file with human protein sequences.

<br />
```{r eval=FALSE}
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
awk '/CHR_HSCHR/{n=2}; n {n--; next}; 1' < human_reference.fa > human_reference_cl.fa
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
<br />
**3. Install Blast**

Install Blast and create a database with the reference sequences.
<br />
```{r eval=FALSE}
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
```{bash}
cd $HOME/conserved_gene_order/reference_human
grep '>' human_reference_cl.fa> ids.txt

```
<br />

Now I have all trnascript ids with their location and the blast results.

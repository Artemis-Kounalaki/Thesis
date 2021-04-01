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
<br /> Find the overlapped proteins (isoforms) that are not mentioned as such and clean the genome from those. <br />
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


# Find overlapped ids and save them

grep '>' human_reference_clean.fa> ids.txt
cd $HOME
python3 ov_ids_new_h.py


# Clean reference genome from overlapped ids

python3 clear_genome_ovrl_h.py
sed -i '/^[[:blank:]]*$/d' ~/conserved_gene_order1/human_reference/human_reference_clean1.fa


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
makeblastdb -in $HOME/conserved_gene_order/human_reference/human_reference_clean1.fa -dbtype prot -parse_seqids  -out database_human


```
**3. Run Blast and find human paralogs.**

Run blastp.
<br />

```
cd $HOME/conserved_gene_order/human_reference
blastp -num_threads 16 -db $HOME/conserved_gene_order1/blast_db/database_human -evalue 1e-6 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen" -qcov_hsp_perc 50 -query $HOME/conserved_gene_order1/human_reference/human_reference_clean1.fa >'results_human.txt'

```

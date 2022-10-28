---
title: "Gene order conservation-sequel- (1)"
author: "Artemis Kounalaki"
---

**Now, its time to find ortholgs.** <br />
## 1st group: Homo sapiens - Macaca mulatta <br />
**1. Download macaca's protein whole genome sequences.**
<br />

```

mkdir $HOME/conserved_gene_order1
cd $HOME/conserved_gene_order1
mkdir human_reference
cd human_reference
wget 'ftp://ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz'
gunzip Homo_sapiens.GRCh38.pep.all.fa.gz


# Make sequences in one line

awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < Homo_sapiens.GRCh38.pep.all.fa > human_reference.fa

```


**2. Clean macaca's protein whole genome sequences.**

It's appropriate to clean file from sequences like scaffold regions, isoforms and unknown regions.<br/>

```
# scaffolds

awk '/scaffold/{n=2}; n {n--; next}; 1' < human_reference.fa > human_reference_cl.fa

# isoforms

awk '/isoform/{n=2}; n {n--; next}; 1' < human_reference_cl.fa > human_reference_cln.fa

# haplotypic regions

awk '/CHR_/{n=2}; n {n--; next}; 1' < human_reference_cln.fa > human_reference_cleann.fa

# MT chromosome

awk '/GRCh38:MT/{n=2}; n {n--; next}; 1' < human_reference_cleann.fa > human_reference_clean.fa

# remove useless files

rm human_reference.fa
rm human_reference_cl.fa
rm human_reference_cln.fa
rm human_reference_cleann.fa


```


**3. Find overlapped ids and save them & clean reference genome from overlapped ids .**
<br/>

```

grep '>' human_reference_clean.fa> ids.txt
cd $HOME
python3 ov_ids_new_h.py
python3 clear_genome_ovrl_h.py
sed -i '/^[[:blank:]]*$/d' ~/conserved_gene_order1/human_reference/human_reference_clean1.fa


```

**4. Make Blast database and run blast.**

Run blast reciprocal to find macaca's proteins that are close enough to human proteins.
<br/>
```
cd $HOME/conserved_gene_order1
#wget https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi-blast-2.10.1+-x64-macosx.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi-blast-2.11.0+-x64-linux.tar.gz
tar -xvzf ncbi-blast-2.11.0+-x64-linux.tar.gz
##tar -xvzf ncbi-blast-2.10.1+-x64-macosx.tar.gz
rm ncbi-blast-2.11.0+-x64-linux.tar.gz
cd $HOME/conserved_gene_order1/ncbi-blast-2.11.0+/bin
export PATH=$PATH:$pwd


# Database

cd ..
mkdir $HOME/conserved_gene_order1/blast_db
export BLASTDB=$HOME/conserved_gene_order1/blast_db
set BLASTDB=$HOME/conserved_gene_order1/blast_db
cd $HOME/conserved_gene_order1/blast_db


# Make the blast database.

makeblastdb -in $HOME/conserved_gene_order1/human_reference/human_reference_clean1.fa -dbtype prot -parse_seqids -out database_human
cd $HOME/conserved_gene_order1/human_reference
# Run blastp.

blastp -num_threads 16 -db $HOME/conserved_gene_order1/blast_db/database_human -evalue 1e-6 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen" -qcov_hsp_perc 50 -query $HOME/conserved_gene_order1/human_reference/human_reference_clean1.fa >'results_human.txt'
cd $HOME/conserved_gene_order1/human_reference

# Now, its time to find ortholgs. 1st group: H.sapiens-macaque

mkdir $HOME/conserved_gene_order1/macaca_reference
cd $HOME/conserved_gene_order1/macaca_reference
wget 'ftp://ftp.ensembl.org/pub/release-100/fasta/macaca_mulatta/pep/Macaca_mulatta.Mmul_10.pep.all.fa.gz'
gunzip Macaca_mulatta.Mmul_10.pep.all.fa.gz

# Make sequences in one line
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < Macaca_mulatta.Mmul_10.pep.all.fa > macaca_reference.fa

# scaffolds

awk '/scaffold/{n=2}; n {n--; next}; 1' < macaca_reference.fa > macaca_reference_cl.fa

# isoforms

awk '/isoform/{n=2}; n {n--; next}; 1' < macaca_reference_cl.fa > macaca_reference_cln.fa

# unkonown regions

awk '/QNVO/{n=2}; n {n--; next}; 1' < macaca_reference_cln.fa > macaca_reference_clean1.fa

# regions with L1 position

awk '/Mmul_10:ML1/{n=2}; n {n--; next}; 1' < macaca_reference_clean1.fa > macaca_reference_clean.fa

# MT chromosome

awk '/Mmul_10:MT/{n=2}; n {n--; next}; 1' < macaca_reference_cleann.fa > macaca_reference_clean.fa

# remove useless files

rm macaca_reference.fa
rm macaca_reference_cl.fa
rm macaca_reference_cln.fa
rm macaca_reference_clean1.fa
rm macaca_reference_cleann.fa


# Find overlapped ids and save them

grep '>' macaca_reference_clean.fa> macaca_ids.txt
cd $HOME
python3 ov_ids_new_m.py

# Clean reference genome from overlapped ids

python3 clear_genome_ovrl_m.py
sed -i '/^[[:blank:]]*$/d' ~/conserved_gene_order1/macaca_reference/macaca_reference_clean1.py

# Make the blast database.

cd $HOME/conserved_gene_order1/blast_db
makeblastdb -in $HOME/conserved_gene_order1/macaca_reference/macaca_reference_clean1.fa -dbtype prot -parse_seqids -out database_macaca
cd $HOME/conserved_gene_order1/macaca_reference


# Run blastp. Human - Macaca

blastp -num_threads 16 -db $HOME/conserved_gene_order1/blast_db/database_macaca -evalue 1e-6 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen" -qcov_hsp_perc 50 -max_hsps 1 -query $HOME/conserved_gene_order1/human_reference/human_reference_clean1.fa >'results_human-macaca.txt'

# Run blastp. Macaca - Human

blastp -num_threads 16 -db $HOME/conserved_gene_order1/blast_db/database_human -evalue 1e-6 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen" -qcov_hsp_perc 50 -max_hsps 1 -query $HOME/conserved_gene_order1/macaca_reference/macaca_reference_clean1.fa >'results_macaca-human.txt'


```


**5. Find reciprocal proteins(human-macaca & macaca-human).**
It's appropriate to exclude proteins that not match reciprocally.

```
# Find reciprocal blast results

cd $HOME
python3 reciprocal_mac-hum.py
cd $HOME/conserved_gene_order1/macaca_reference
cat reciprocal_hum-mac*.txt >> reciprocal_hum-mac_all.txt
tr -d "['],"< reciprocal_hum-mac_all.txt > reciprocal_h-m.txt
rm reciprocal_hum-mac*.txt
 ```
<br/>
 **6. Find best hits.**
```
cd $HOME
python3 mybl_1.py
cd $HOME
```
<br/>
**7. Find GCOs and statistics.**
```
python3 order_hum-mac.py
python3 CGO_res_h-m.py
python3 statistics_h-mac.py
```

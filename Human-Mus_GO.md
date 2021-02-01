---
title: "Gene order conservation-sequel- (2)"
author: "Artemis Kounalaki"
---

## 2d group: Homo sapiens - Mus musculus <br />
The same pipeline as Human - Macaca was used for compare Human - Mus musculus
*1. Download macaca's protein whole genome sequences.**<br />
```
mkdir $HOME/conserved_gene_order/mus_reference
cd $HOME/conserved_gene_order/mus_reference
wget 'ftp://ftp.ensembl.org/pub/release-100/fasta/mus_musculus/pep/Mus_musculus.GRCm38.pep.all.fa.gz'
gunzip Mus_musculus.GRCm38.pep.all.fa.gz


# Make its sequence in one line

awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < Mus_musculus.GRCm38.pep.all.fa > mus_reference.fa

```
**2. Clean mus' protein whole genome sequences.**

It's appropriate to clean file from sequences like scaffold regions, isoforms and unknown regions.<br/>

```
# scaffolds

awk '/scaffold/{n=2}; n {n--; next}; 1' < mus_reference.fa > mus_reference_cl.fa

# isoforms

awk '/isoform/{n=2}; n {n--; next}; 1' < mus_reference_cl.fa > mus_reference_clean.fa

# Remove useless files

rm mus_reference.fa
rm mus_reference_cl.fa
```

**3. Extract protein ids and their location.**
<br/>
```
# Find ids

grep '>' mus_reference_clean.fa> mus_ids.txt
```
**4. Make Blast database and run blast.**

Run blast reciprocal to find mus' proteins that are close enough to human proteins.
<br/>
```
# Make the blast database.

cd $HOME/conserved_gene_order/blast_db
makeblastdb -in $HOME/conserved_gene_order/mus_reference/mus_reference_clean.fa -dbtype prot -parse_seqids -out database_mus
cd $HOME/conserved_gene_order/mus_reference
# Run blastp. Human - Mus_musculus

blastp -num_threads 16 -db $HOME/conserved_gene_order/blast_db/database_mus -evalue 1e-10 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen" -qcov_hsp_perc 80 -query $HOME/conserved_gene_order/human_reference/human_reference_clean.fa >'results_human-mus.txt'


# Run blastp. Mus_musculus - Human

blastp -num_threads 16 -db $HOME/conserved_gene_order/blast_db/database_human -evalue 1e-10 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen" -qcov_hsp_perc 80 -query $HOME/conserved_gene_order/mus_reference/mus_reference_clean.fa >'results_mus-human.txt'


# Run blastp. Mus_musculus - Mus_musculus

blastp -num_threads 16 -db $HOME/conserved_gene_order/blast_db/database_mus -evalue 1e-10 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen" -qcov_hsp_perc 80 -query $HOME/conserved_gene_order/mus_reference/mus_reference_clean.fa >'results_mus.txt'

cd $HOME
```
**4. Clean mus blastp results from overlapped ids, clean reciporcal results from overlapped ids and find CGO and nCGO.**<br />
```
python3 overlapped_ids_mus.py
python3 clean_ov_ids_mus.py


# Find reciprocal blast results

python3 reciprocal_mus-hum.py
cd $HOME/conserved_gene_order/mus_reference
cat reciprocal_hum-mus*.txt >> reciprocal_hum-mus_all.txt
tr -d "['],"< reciprocal_hum-mus_all.txt > reciprocal_h-mus.txt
rm reciprocal_hum-mus*.txt
cd $HOME
python3 clean_ov_ids_h-mus.py
python3 order2.py
```

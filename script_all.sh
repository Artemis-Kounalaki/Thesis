# Make a file for the whole process.

mkdir $HOME/conserved_gene_order
cd $HOME/conserved_gene_order


# Make a file inside with the human protein sequences from Ensembl.

mkdir human_reference
cd human_reference
wget 'ftp://ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz'
gunzip Homo_sapiens.GRCh38.pep.all.fa.gz
head -2000 Homo_sapiens.GRCh38.pep.all.fa > Homo_sapiens.GRCh38.pep.all1.fa

# Make sequences in one line

awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < Homo_sapiens.GRCh38.pep.all1.fa > human_reference.fa

# scaffolds

awk '/scaffold/{n=2}; n {n--; next}; 1' < human_reference.fa > human_reference_cl.fa

# isoforms

awk '/isoform/{n=2}; n {n--; next}; 1' < human_reference_cl.fa > human_reference_cln.fa

# haplotypic regions

awk '/CHR_/{n=2}; n {n--; next}; 1' < human_reference_cln.fa > human_reference_clean.fa

# remove useless files

rm human_reference.fa
rm human_reference_cl.fa
rm human_reference_cln.fa



# Blast installation from ncbi

cd $HOME/conserved_gene_order
#wget https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi-blast-2.10.1+-x64-macosx.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi-blast-2.11.0+-x64-linux.tar.gz
tar -xvzf ncbi-blast-2.11.0+-x64-linux.tar.gz
#tar -xvzf ncbi-blast-2.10.1+-x64-macosx.tar.gz
rm ncbi-blast-2.11.0+-x64-linux.tar.gz
cd $HOME/conserved_gene_order/ncbi-blast-2.11.0+/bin
export PATH=$PATH:$pwd


# Database

cd ..
mkdir $HOME/conserved_gene_order/blast_db
export BLASTDB=$HOME/conserved_gene_order/blast_db
set BLASTDB=$HOME/conserved_gene_order/blast_db
cd $HOME/conserved_gene_order/blast_db


# Make the blast database.

makeblastdb -in $HOME/conserved_gene_order/human_reference/human_reference_clean.fa -dbtype prot -parse_seqids -out database_human
cd $HOME/conserved_gene_order/human_reference


# Run blastp.

blastp -num_threads 16 -db $HOME/conserved_gene_order/blast_db/database_human -evalue 1e-10 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen" -qcov_hsp_perc 80 -query $HOME/conserved_gene_order/human_reference/human_reference_clean.fa >'results_human.txt'
cd $HOME/conserved_gene_order/human_reference


# Find the ids

grep '>' human_reference_clean.fa> ids.txt
cd $HOME
python3 overlapped_ids_h.py
python3 clean_ov_ids_h.py





# Now, its time to find ortholgs. 1st group: H.sapiens-macaque

mkdir $HOME/conserved_gene_order/macaca_reference
cd $HOME/conserved_gene_order/macaca_reference
wget 'ftp://ftp.ensembl.org/pub/release-100/fasta/macaca_mulatta/pep/Macaca_mulatta.Mmul_10.pep.all.fa.gz'
gunzip Macaca_mulatta.Mmul_10.pep.all.fa.gz
head -2000 Macaca_mulatta.Mmul_10.pep.all.fa > Macaca_mulatta.Mmul_10.pep.all1.fa

# Make sequences in one line
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < Macaca_mulatta.Mmul_10.pep.all1.fa > macaca_reference.fa

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

# Find ids

grep '>' macaca_reference_clean.fa> macaca_ids.txt

# Make the blast database.

cd $HOME/conserved_gene_order/blast_db
makeblastdb -in $HOME/conserved_gene_order/macaca_reference/macaca_reference_clean.fa -dbtype prot -parse_seqids -out database_macaca
cd $HOME/conserved_gene_order/macaca_reference


# Run blastp. Human - Macaca

blastp -num_threads 16 -db $HOME/conserved_gene_order/blast_db/database_macaca -evalue 1e-10 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen" -qcov_hsp_perc 80 -query $HOME/conserved_gene_order/human_reference/human_reference_clean.fa >'results_human-macaca.txt'


# Run blastp. Macaca - Human

blastp -num_threads 16 -db $HOME/conserved_gene_order/blast_db/database_human -evalue 1e-10 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen" -qcov_hsp_perc 80 -query $HOME/conserved_gene_order/macaca_reference/macaca_reference_clean.fa >'results_macaca-human.txt'


# Run blastp. Macaca-Macaca

blastp -num_threads 16 -db $HOME/conserved_gene_order/blast_db/database_macaca -evalue 1e-10 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen" -qcov_hsp_perc 80 -query $HOME/conserved_gene_order/macaca_reference/macaca_reference_clean.fa >'results_macaca.txt'


# Clean macaca blast results from isoforms.

cd $HOME
python3 overlapped_ids_m.py
python3 clean_ov_ids_m.py

# Find reciprocal blast results

python3 reciprocal_mac-hum.py
cd $HOME/conserved_gene_order/macaca_reference
cat reciprocal_hum-mac*.txt >> reciprocal_hum-mac_all.txt
tr -d "['],"< reciprocal_hum-mac_all.txt > reciprocal_h-m.txt
rm reciprocal_hum-mac*.txt
cd $HOME
python3 clean_ov_ids_h-mac.py
python3 order1.py


# 2st group: H.sapiens- Mus musculus

mkdir $HOME/conserved_gene_order/mus_reference
cd $HOME/conserved_gene_order/mus_reference
wget 'ftp://ftp.ensembl.org/pub/release-100/fasta/mus_musculus/pep/Mus_musculus.GRCm38.pep.all.fa.gz'
gunzip Mus_musculus.GRCm38.pep.all.fa.gz
head -2000 Mus_musculus.GRCm38.pep.all.fa > Mus_musculus.GRCm38.pep.all1.fa
# Make its sequence in one line

awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < Mus_musculus.GRCm38.pep.all1.fa > mus_reference.fa

# scaffolds

awk '/scaffold/{n=2}; n {n--; next}; 1' < mus_reference.fa > mus_reference_cl.fa

# isoforms

awk '/isoform/{n=2}; n {n--; next}; 1' < mus_reference_cl.fa > mus_reference_clean.fa

# Remove useless files

rm mus_reference.fa
rm mus_reference_cl.fa

# Find ids

grep '>' mus_reference_clean.fa> mus_ids.txt

# Make the blast database.

cd $HOME/conserved_gene_order/blast_db
makeblastdb -in $HOME/conserved_gene_order/mus_reference/mus_reference_clean.fa -dbtype prot -parse_seqids -out database_mus
cd $HOME/conserved_gene_order/mus_reference


# Run blastp. Mus musculus - Human

blastp -num_threads 16 -db $HOME/conserved_gene_order/blast_db/database_human -evalue 1e-10 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen" -qc$

# Run bastp. Human - Mus musculus

blastp -num_threads 16 -db $HOME/conserved_gene_order/blast_db/database_mus -evalue 1e-10 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen" -qcov$

# Run blastp. Mus - Mus

blastp -num_threads 16 -db $HOME/conserved_gene_order/blast_db/database_mus -evalue 1e-10 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen" -qcov$

cd $HOME
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

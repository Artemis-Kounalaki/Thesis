# Make a file for the whole process.

mkdir $HOME/conserved_gene_order1
cd $HOME/conserved_gene_order1


# Make a file inside with the human protein sequences from Ensembl.

mkdir human_reference
cd human_reference
wget 'ftp://ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz'
gunzip Homo_sapiens.GRCh38.pep.all.fa.gz


# Make sequences in one line

awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < Homo_sapiens.GRCh38.pep.all.fa > human_reference.fa

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

# Find overlapped ids and save them

grep '>' human_reference_clean.fa> ids.txt
cd $HOME
python3 ov_ids_new_h.py

# Clean reference genome from overlapped ids

python3 clear_genome_ovrl_h.py
sed -i '/^[[:blank:]]*$/d' ~/conserved_gene_order1/human_reference/human_reference_clean1.fa


# Blast installation from ncbi

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


# Run blastp. Macaca-Macaca

#blastp -num_threads 16 -db $HOME/conserved_gene_order/blast_db/database_macaca -evalue 1e-10 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen" -qcov_hsp_perc 80 -query $HOME/conserved_gene_order/macaca_reference/macaca_reference_clean.fa >'results_macaca.txt'


# Find reciprocal blast results

cd $HOME
python3 reciprocal_mac-hum.py
cd $HOME/conserved_gene_order1/macaca_reference
cat reciprocal_hum-mac*.txt >> reciprocal_hum-mac_all.txt
tr -d "['],"< reciprocal_hum-mac_all.txt > reciprocal_h-m.txt
rm reciprocal_hum-mac*.txt
cd $HOME
python3 mybl_1.py
cd $HOME
# ena stop edw
python3 order1.py
cd $HOME/conserved_gene_order1/macaca_reference
python3 CGO_res.py

# 2st group: H.sapiens- Mus musculus

mkdir $HOME/conserved_gene_order1/mus_reference
cd $HOME/conserved_gene_order1/mus_reference
wget 'ftp://ftp.ensembl.org/pub/release-100/fasta/mus_musculus/pep/Mus_musculus.GRCm38.pep.all.fa.gz'
gunzip Mus_musculus.GRCm38.pep.all.fa.gz

# Make its sequence in one line

awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < Mus_musculus.GRCm38.pep.all.fa > mus_reference.fa

# scaffolds

awk '/scaffold/{n=2}; n {n--; next}; 1' < mus_reference.fa > mus_reference_cl.fa

# isoforms

awk '/isoform/{n=2}; n {n--; next}; 1' < mus_reference_cl.fa > mus_reference_clean.fa

# Remove useless files

rm mus_reference.fa
rm mus_reference_cl.fa

# Find overlapped ids and save them

grep '>' mus_reference_clean.fa> mus_ids.txt
cd $HOME
python3 ov_ids_new_mus.py

# Clean reference genome from overlapped ids

python3 clear_genome_ovrl_mus.py
sed -i '/^[[:blank:]]*$/d' ~/conserved_gene_order1/mus_reference/mus_reference_clean1.fa


# Make the blast database.

cd $HOME/conserved_gene_order1/blast_db
makeblastdb -in $HOME/conserved_gene_order1/mus_reference/mus_reference_clean1.fa -dbtype prot -parse_seqids -out database_mus
cd $HOME/conserved_gene_order1/mus_reference


# Run blastp. Mus musculus - Human

blastp -num_threads 16 -db $HOME/conserved_gene_order1/blast_db/database_human -evalue 1e-6 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen" -qcov_hsp_perc 50 -max_hsps 1 -max_target_seqs 1 -query $HOME/conserved_gene_order1/mus_reference/mus_reference_clean1.fa >'results_mus-human.txt'

# Run bastp. Human - Mus musculus

blastp -num_threads 16 -db $HOME/conserved_gene_order1/blast_db/database_mus -evalue 1e-6 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen"  -qcov_hsp_perc 50 -max_hsps 1 -max_target_seqs 1 -query $HOME/conserved_gene_order1/human_reference/human_reference_clean1.fa >'results_human-mus.txt'

# Run blastp. Mus - Mus

#blastp -num_threads 16 -db $HOME/conserved_gene_order/blast_db/database_mus -evalue 1e-10 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen" -qcov_hsp_perc 80 -max_hsps 1 -max_target_seqs 1 -query $HOME/conserved_gene_order/mus_reference/mus_reference_clean.fa >'results_mus.txt'


# Find reciprocal blast results

cd $HOME
python3 reciprocal_mus-hum.py
cd $HOME/conserved_gene_order1/mus_reference
cat reciprocal_hum-mus*.txt >> reciprocal_hum-mus_all.txt
tr -d "['],"< reciprocal_hum-mus_all.txt > reciprocal_h-mus.txt
rm reciprocal_hum-mus*.txt
cd $HOME
python3 order2.py

# 3d Group : Mus- Macaca

cd $HOME/conserved_gene_order1/mus_reference

# Run blastp. Mus - Macaca

blastp -num_threads 16 -db $HOME/conserved_gene_order1/blast_db/database_mus -evalue 1e-6 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen" -qcov_$


# Run blastp. Macaca - Mus

blastp -num_threads 16 -db $HOME/conserved_gene_order1/blast_db/database_macaca -evalue 1e-6 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen" -qc$


# Find reciprocal blast results

cd $HOME
python3 reciprocal_mac-mus.py
cd $HOME/conserved_gene_order1/mus_reference
cat reciprocal_mus-mac*.txt >> reciprocal_mus-mac_all.txt
tr -d "['],"< reciprocal_mus-mac_all.txt > reciprocal_mus-m.txt
rm reciprocal_mus-mac*.txt
cd $HOME
python3 order3.py

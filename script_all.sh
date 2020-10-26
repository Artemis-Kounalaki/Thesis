# Make a file for the whole process.

mkdir $HOME/conserved_gene_order
cd $HOME/conserved_gene_order


# Make a file inside with the human protein sequences from Ensembl.

mkdir human_reference
cd human_reference
wget 'ftp://ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz'
gunzip Homo_sapiens.GRCh38.pep.all.fa.gz


#make sequences in one line

awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < Homo_sapiens.GRCh38.pep.all.fa > human_reference.fa


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

makeblastdb -in $HOME/conserved_gene_order/human_reference/human_reference_clean.fa -dbtype prot -parse_seqids -out database_human
cd $HOME/conserved_gene_order/human_reference


# Run blastp.

blastp -num_threads 16 -db $HOME/conserved_gene_order/blast_db/database_human -evalue 1e-10 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen" -qcov_hsp_perc 80 -query $HOME/conserved_gene_order/human_reference/human_reference_clean.fa >'results_human.txt'
cd $HOME/conserved_gene_order/human_reference


# Find the ids

grep '>' human_reference_clean.fa> ids.txt
cd $HOME
python3 overlapped_ids.py
python3 clean_ovrl.py

cd $HOME/conserved_gene_order/human_reference
tr -d "['],"< clean_blast_res.txt > clean_blast_results.txt
rm clean_blast_res.txt





# Now, its time to find ortholgs. 1st group: H.sapiens-macaque

mkdir $HOME/conserved_gene_order/macaca_reference
cd $HOME/conserved_gene_order/macaca_reference
wget 'ftp://ftp.ensembl.org/pub/release-100/fasta/macaca_mulatta/pep/Macaca_mulatta.Mmul_10.pep.all.fa.gz'
gunzip Macaca_mulatta.Mmul_10.pep.all.fa.gz

#make sequences in one line
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < Macaca_mulatta.Mmul_10.pep.all.fa > macaca_reference.fa

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


#Find ids
grep '>' macaca_reference_clean.fa> macaca_ids.txt

# Make the blast database.

cd $HOME/conserved_gene_order/blast_db
makeblastdb -in $HOME/conserved_gene_order/macaca_reference/macaca_reference_clean.fa -dbtype prot -parse_seqids -out database_macaca
cd $HOME/conserved_gene_order/macaca_reference


# Run blastp. Human - Macaca

blastp -num_threads 16 -db $HOME/conserved_gene_order/blast_db/database_macaca -evalue 1e-10 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen" -qcov_hsp_perc 80 -query $HOME/conserved_gene_order/human_reference/human_reference_clean.fa >'results_human-macaca.txt'


# Run blastp. Macaca - Human

blastp -num_threads 16 -db $HOME/conserved_gene_order/blast_db/database_human -evalue 1e-10 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen" -qcov_hsp_perc 80 -query $HOME/conserved_gene_order/macaca_reference/macaca_reference_clean.fa >'results_macaca-human.txt'

cd $HOME
python3 reciprocal_mac-hum.py
cd $HOME/conserved_gene_order/macaca_reference
cat reciprocal_hum-mac*.txt >> reciprocal_hum-mac_all.txt
tr -d "['],"< reciprocal_hum-mac_all.txt > reciprocal_h-m.txt
rm reciprocal_hum-mac.*.txt


#  Make groups: human proteins that are similar with macaca's proteins.

cd $HOME
python3 protein_groups.py

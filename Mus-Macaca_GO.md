---
title: "Gene order conservation-sequel- (1)"
author: "Artemis Kounalaki"
---


## 3d group: Macaca mulatta - Mus musculus <br />
**1. Run blastp between the two organisms.**
<br />

cd $HOME/conserved_gene_order1/mus_reference

# Run blastp. Mus - Macaca

blastp -num_threads 16 -db $HOME/conserved_gene_order1/blast_db/database_mus -evalue 1e-6 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen" -qcov_hsp_perc 50 -max_hsps 1 -max_target_seqs 1 -query $HOME/conserved_gene_order1/macaca_reference/macaca_reference_clean1.fa > 'results_mus-macaca.txt'


# Run blastp. Macaca - Mus

blastp -num_threads 16 -db $HOME/conserved_gene_order1/blast_db/database_macaca -evalue 1e-6 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen" -qcov_hsp_perc 50 -max_hsps 1 -max_target_seqs 1 -query $HOME/conserved_gene_order1/mus_reference/mus_reference_clean1.fa >'results_macaca-mus.txt'

**1. Find reciprocal proteins & conserved -non conserved gene order proteins.**
# Find reciprocal blast results

cd $HOME
python3 reciprocal_mac-mus.py
cd $HOME/conserved_gene_order1/mus_reference
cat reciprocal_mus-mac*.txt >> reciprocal_mus-mac_all.txt
tr -d "['],"< reciprocal_mus-mac_all.txt > reciprocal_mus-m.txt
rm reciprocal_mus-mac*.txt
cd $HOME
python3 order3.py

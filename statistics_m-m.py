import pandas as pd
from scipy import stats
import os
import warnings

# Function with statistics in GCOS/nGCOs in group 3: mouse-macaca

def stat(path,txt_file, reciprocal_file):

    # Statistic between GCOs and nGCOs identity
    os.chdir(os.path.expanduser(path))
    warnings.simplefilter(action='ignore', category=FutureWarning)
    pd.options.mode.chained_assignment = None

    f= open(txt_file,"w+")
    f.write('This is a statistical test applying on the results of RBH between groups of GCOs/nGCOs.''\n')

    # Load GCOs nGCOs results
    df_al = pd.read_csv(reciprocal_file, sep='\t')
    df_al.columns=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'slen', 'qlen', 'Conserved']
    identity_cgo = df_al.pident[df_al.Conserved==1]
    identity_ncgo = df_al.pident[df_al.Conserved==0]
    f.write('Two-sample Mann–Whitney U test applied between percentage identity of GCOs and nGCOs...''\n')
    f.write('The number of GCOs is: %f \n' %len(identity_cgo))
    f.write('The number of nGCOs is : %f \n'  %len(identity_ncgo))
    f.write('The average perc identity of GCOs is: %f \n'% (sum(identity_cgo)/len(identity_cgo)))
    f.write('The average perc identity of nGCOs is : %f \n '% (sum(identity_ncgo)/len(identity_ncgo)))
    f.write('T-statistic : %d and p-value : %g in GCOs & nGCOs comparison of percentage identity. \n' %(stats.mannwhitneyu(identity_cgo, identity_ncgo)))

stat('~/conserved_gene_order1/mus_reference', 'statistics_m-m.txt' ,'CGO_nCGO_m-m.txt')

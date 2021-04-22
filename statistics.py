import pandas as pd
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def statistics(path1, file, path_results_paralogs,file_results_paralogs, save_name_cgo,save_name_ncgo):
    # Statistic between cgo and ncgo identity
    os.chdir(os.path.expanduser(path1))
    print('This is a statistical test applying on the results of bastp between groups and is about hits.')
    df_al = pd.read_csv(file, sep='\t')
    identity_cgo = df_al.P_identity[df_al.Conserved==1]
    identity_ncgo = df_al.P_identity[df_al.Conserved==0]
    print('T-test between identity percentage of cgo and ncgo...''\n')
    print('The number of cgo genes is:', len(identity_cgo))
    print('The number of non cgo genes is:', len(identity_ncgo))
    print('The average identity of cgo is:',sum(identity_cgo)/len(identity_cgo))
    print('The average identity of non cgo is:',sum(identity_ncgo)/len(identity_ncgo))
    print(stats.ttest_ind(identity_cgo, identity_ncgo, equal_var = False))

    # Find paralogs
    print('Loading human paralogs ...')

    os.chdir(os.path.expanduser(path_results_paralogs))
    df_paralogs = pd.read_csv(file_results_paralogs, sep='\t', header=None)
    df_paralogs.columns=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'slen', 'qlen']
    df_paralogs.qseqid = df_paralogs.qseqid .astype(str)
    df_paralogs.sseqid = df_paralogs.sseqid .astype(str)
    df_del = df_paralogs[df_paralogs.qseqid == df_paralogs.sseqid] # delete the hit of a protein with itself
    df_paralogs=df_paralogs[~df_paralogs.isin(df_del)]
    df_paralogs = df_paralogs.dropna(how='all')
    df_paralogs = df_paralogs.reset_index(drop=True)
    df_paralogs["pairs"] = df_paralogs.apply(lambda row: ''.join(sorted([row["qseqid"], row["sseqid"]])), axis = 1)
    only_pairs = df_paralogs[df_paralogs["pairs"].duplicated(keep = False)].sort_values(by = "pairs") # find reciprocallity
    only_pairs = only_pairs.reset_index(drop=True)


    # Statistic paralog conserved and paralog non conserved

    print('Statistical analysis between subgroups is loading...')
    df_all = pd.read_csv(file, sep='\t')
    df_all.columns=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'slen', 'qlen', 'Conserved']
    genes=df_all[['qseqid','sseqid']].values.ravel()
    par= pd.unique(only_pairs[['qseqid','sseqid']].values.ravel())

    print('The number of human paralogs is',len(par))
    df_all['Paralog']= df_all.qseqid.isin(par)
    df_all['Paralog']= df_all['Paralog'].fillna(0)
    df_all['Paralog']= df_all['Paralog'].astype(int)
    cons= df_all[df_all.Conserved==1]
    cons.pident = cons.pident.astype(float)
    ncons= df_all[df_all.Conserved==0]
    ncons.pident = ncons.pident.astype(float)

    print('The number of conserved genes that are paralogs is:')
    print(len(cons[cons.qseqid.isin(par)]))
    print('Their average percentage identity is:',sum(cons.pident[cons.qseqid.isin(par)])/len(cons.pident[cons.qseqid.isin(par)]))

    print('The number of non conserved genes that are paralogs is:')
    print(len(ncons[ncons.qseqid.isin(par)]))
    print('Their average percentage isentity is:',sum(ncons.pident[ncons.qseqid.isin(par)])/len(ncons.pident[ncons.qseqid.isin(par)]))
    print('T- test between cgo paralogs & ncgo paralogs')
    print(stats.ttest_ind(cons.pident[cons.qseqid.isin(par)], ncons.pident[ncons.qseqid.isin(par)], equal_var = False))

    cons1=cons[~cons.qseqid.isin(par)]
    print('The number of conserved genes without paralogs is:',len(cons1))
    ncons1 = ncons[~ncons.qseqid.isin(par)]
    print('The number of non conserved genes without paralogs is:',len(ncons1))

    print('T-test between percentage identity of non conserved genes that are paralogs and those that are not paralogs...')
    print(stats.ttest_ind(ncons.pident[ncons.qseqid.isin(par)], ncons1.pident, equal_var = False))
    print('T-test between percentage identity of conserved genes that are paralogs and those that are not paralogs...')
    print(stats.ttest_ind(cons.pident[cons.qseqid.isin(par)], cons1.pident, equal_var = False))
    print('t-test between cgo non paralogs and ncgo non paralogs')
    print(stats.ttest_ind(ncons1.pident, cons1.pident, equal_var = False))

    cons.qseqid[cons.qseqid.isin(par)].to_csv(save_name_cgo, sep='\t',index=None)
    ncons.qseqid[ncons.qseqid.isin(par)].to_csv(save_name_ncgo, sep='\t',index=None)



    pal = sns.color_palette('Paired')
    sns.boxplot(x='Paralog', y='pident', hue='Paralog', data=cons,   palette='husl')

    sns.boxplot(x='Conserved', y='pident', hue='Conserved', data=df_all,   palette='husl')

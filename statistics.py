import pandas as pd
from scipy import stats
import seaborn as sns
from matplotlib import pyplot
import warnings

def stat(path1,txt_file, reciprocal_file, path2):

    # Statistic between GCOs and nGCOs identity
    os.chdir(os.path.expanduser(path1))
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


    # Find paralogs
    f.write('Loading human paralogs ... \n')

    os.chdir(os.path.expanduser(path2))
    df_paralogs = pd.read_csv('results_human.txt', sep='\t', header=None)
    df_paralogs.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
                           'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'slen', 'qlen']
    df_paralogs.qseqid = df_paralogs.qseqid .astype(str)
    df_paralogs.sseqid = df_paralogs.sseqid .astype(str)


    # Delete the hit of a protein with itself

    df_del = df_paralogs[df_paralogs.qseqid == df_paralogs.sseqid]
    df_paralogs = df_paralogs[~df_paralogs.isin(df_del)]
    df_paralogs = df_paralogs.dropna(how='all')
    df_paralogs = df_paralogs.reset_index(drop=True)


    # Reciprocals
    df_paralogs["pairs"] = df_paralogs.apply(lambda row: ''.join(sorted([row["qseqid"], row["sseqid"]])), axis=1)
    only_pairs = df_paralogs[df_paralogs["pairs"].duplicated(keep = False)].sort_values(by = "pairs") # find reciprocallity
    keep=list(set(only_pairs['qseqid']).intersection(set(only_pairs['sseqid']))) #prevent false pairs ,alphabetically sort pairs
    only_pairs=only_pairs[only_pairs['qseqid'].isin(keep) & only_pairs['sseqid'].isin(keep)]
    only_pairs = only_pairs.reset_index(drop=True)
    par = pd.unique(only_pairs[['qseqid', 'sseqid']].values.ravel())
    os.chdir(os.path.expanduser(path1))


    # Statistic paralog conserved and paralog non conserved
    f.write('Statistical analysis between subgroups is loading... \n')
    f.write('The number of human paralogs is: %d \n '%len(par))
    df_al['Paralog']= df_al.qseqid.isin(par)
    df_al['Paralog']= df_al['Paralog'].fillna(0)
    df_al['Paralog']= df_al['Paralog'].astype(int)
    cons= df_al[df_al.Conserved==1]
    cons.pident = cons.pident.astype(float)
    ncons= df_al[df_al.Conserved==0]
    ncons.pident = ncons.pident.astype(float)

    f.write('The number of GCOs that have a paralog is:')
    f.write('%f \n' %len(cons[cons.qseqid.isin(par)]))
    f.write('Their average percentage identity is: %f \n'%(sum(cons.pident[cons.qseqid.isin(par)])/len(cons.pident[cons.qseqid.isin(par)])))
    f.write('The number of nGCOs that have a paralog is:')
    f.write('%f \n' %len(ncons[ncons.qseqid.isin(par)]))
    f.write('Their average percentage identity is: %f \n'%(sum(ncons.pident[ncons.qseqid.isin(par)])/len(ncons.pident[ncons.qseqid.isin(par)])))
    f.write('\n Two sample Mann–Whitney U test between GCOs with paralog & nGCOs with paralog:')
    f.write('\n T-statistic & p-value: %g, %g' %(stats.mannwhitneyu(cons.pident[cons.qseqid.isin(par)], ncons.pident[ncons.qseqid.isin(par)])))

    cons1=cons[~cons.qseqid.isin(par)]
    f.write('\n The number of GCOs without paralogs is: %d'%len(cons1))
    ncons1 = ncons[~ncons.qseqid.isin(par)]
    f.write('\n The number of nGCOs without paralogs is: %d \n'%len(ncons1))
    f.write('The average percentage identity of GCOs without paralog is: %f '%(sum(cons1.pident)/len(cons1)))
    f.write('\n The average percentage identity of nGCOs without paralog is: %f \n'%(sum(ncons1.pident)/len(ncons1)))
    f.write('Two sample Mann–Whitney U test between percentage identity of nGCOs that have a paralog and those not have: ')
    f.write('%g, %g \n' %(stats.mannwhitneyu(ncons.pident[ncons.qseqid.isin(par)], ncons1.pident)))
    f.write('Two sample Mann–Whitney U test between percentage identity of GCOs that have a paralog and those that not have: ')
    f.write('%g, %g \n' %(stats.mannwhitneyu(cons.pident[cons.qseqid.isin(par)], cons1.pident)))
    f.write('Two sample Mann–Whitney U test between GCOs not having paralog and nGCOs not having paralog: ')
    f.write('%g, %g \n' %(stats.mannwhitneyu(ncons1.pident, cons1.pident)))
    f.close()

    #cons.qseqid[cons.qseqid.isin(par)].to_csv('cons_par1.txt', sep='\t',index=None)
    #ncons.qseqid[ncons.qseqid.isin(par)].to_csv('ncons_par1.txt', sep='\t',index=None)

    #Create a full df with 2 columns , pident and target for dist plot

    df_a=pd.DataFrame(cons.pident[~cons.qseqid.isin(par)],columns = ["pident"])
    df_a=df_a.reset_index(drop=True)
    df_a['target']=0
    df_a1=pd.DataFrame(cons.pident[cons.qseqid.isin(par)],columns = ["pident"])
    df_a1=df_a1.reset_index(drop=True)
    df_a1['target']=1
    df_a2=pd.DataFrame(ncons.pident[~ncons.qseqid.isin(par)],columns = ["pident"])
    df_a2=df_a2.reset_index(drop=True)
    df_a2['target']=2
    df_a3=pd.DataFrame(ncons.pident[ncons.qseqid.isin(par)],columns = ["pident"])
    df_a3=df_a3.reset_index(drop=True)
    df_a3['target']=3
    df_aa = pd.concat([df_a, df_a1, df_a2, df_a3])
    df_aa=df_aa.reset_index(drop=True)


    # Distribution plot 4 groups h-m / h-mus

    sns.set_palette(palette="Set2",n_colors=4)
    a4_dims = (10, 8)
    sns.set_style('white')
    fig, ax = pyplot.subplots(figsize=a4_dims)
    sns.distplot(cons.pident[~cons.qseqid.isin(par)], hist=False, rug=True, kde_kws=dict(linewidth=1.25))
    sns.distplot(cons.pident[cons.qseqid.isin(par)], hist=False, rug=True,kde_kws=dict(linewidth=1.25))
    sns.distplot(ncons.pident[~ncons.qseqid.isin(par)], hist=False, rug=True,kde_kws=dict(linewidth=1.25))
    sns.distplot(ncons.pident[ncons.qseqid.isin(par)], hist=False, rug=True,kde_kws=dict(linewidth=1.25))
    ax.set_title('Distribution plot')
    ax.set_xlabel('Percentage Identity')
    fig.legend(labels=["GCOs in single copy","GCOs with paralog","nGCOs in single copy","nGCOs with paralog"])

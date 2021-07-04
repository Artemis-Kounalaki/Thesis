import pandas as pd
from scipy import stats
from scipy.stats import mannwhitneyu
import seaborn as sns
from matplotlib import pyplot


# Statistic between cgo and ncgo identity


#os.chdir(os.path.expanduser(path1))
# Open a txt file and plcae statistical results.
f= open("testanal.txt","w+")

f.write('This is a statistical test applying on the results of bastp between groups and is about hits.''\n')
df_al = pd.read_csv('al_mus.txt', sep='\t')
identity_cgo = df_al.P_identity[df_al.Conserved==1]
identity_ncgo = df_al.P_identity[df_al.Conserved==0]
f.write('Two-sample wilcoxon test between identity percentage of cgo and ncgo...''\n')
f.write('The number of cgo genes is: %f \n' %len(identity_cgo))
f.write('The number of non cgo genes is : %f \n'  %len(identity_ncgo))
f.write('The average identity of cgo is: %f \n'% (sum(identity_cgo)/len(identity_cgo)))
f.write('The average identity of non cgo is : %f \n '%sum((identity_ncgo)/len(identity_ncgo)))
f.write('T-statistic : %d and p-value : %g in cgo & non cgo comparison \n' %(stats.mannwhitneyu(identity_cgo, identity_ncgo)))

# Find paralogs
f.write('Loading human paralogs ... \n')

#os.chdir(os.path.expanduser(path_results_paralogs))
df_paralogs = pd.read_csv('results_human.txt', sep='\t', header=None)
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

f.write('Statistical analysis between subgroups is loading... \n')
df_all = pd.read_csv('al_mus.txt', sep='\t')
df_all.columns=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'slen', 'qlen', 'Conserved']
genes=df_all[['qseqid','sseqid']].values.ravel()
par= pd.unique(only_pairs[['qseqid','sseqid']].values.ravel())

f.write('The number of human paralogs is: %d \n '%len(par))
df_all['Paralog']= df_all.qseqid.isin(par)
df_all['Paralog']= df_all['Paralog'].fillna(0)
df_all['Paralog']= df_all['Paralog'].astype(int)
cons= df_all[df_all.Conserved==1]
cons.pident = cons.pident.astype(float)
ncons= df_all[df_all.Conserved==0]
ncons.pident = ncons.pident.astype(float)

f.write('The number of conserved genes that are paralogs is:')
f.write('%f' %len(cons[cons.qseqid.isin(par)]))
f.write('Their average percentage identity is: %f \n'%(sum(cons.pident[cons.qseqid.isin(par)])/len(cons.pident[cons.qseqid.isin(par)])))

f.write('The number of non conserved genes that are paralogs is:')
f.write('%f \n' %len(ncons[ncons.qseqid.isin(par)]))
f.write('Their average percentage isentity is: %f \n'%(sum(ncons.pident[ncons.qseqid.isin(par)])/len(ncons.pident[ncons.qseqid.isin(par)])))
f.write('Two sample wilcoxon test between CGOs with paralog & non-CGOs with paralog')
f.write('%g, %g \n' %(stats.mannwhitneyu(cons.pident[cons.qseqid.isin(par)], ncons.pident[ncons.qseqid.isin(par)])))

cons1=cons[~cons.qseqid.isin(par)]
f.write('The number of conserved genes without paralogs is: %d'%len(cons1))
ncons1 = ncons[~ncons.qseqid.isin(par)]
f.write('The average percentage identity of conserved genes without paralogs is %f:'%(sum(cons1.pident)/len(cons1)))


f.write('The number of non conserved genes without paralogs is: %d \n'%len(ncons1))
f.write('The average percentage identity of non conserved genes without paralog is:%f \n'%(sum(ncons1.pident)/len(ncons1)))
f.write('Two sample wilcoxon test between percentage identity of non conserved genes that are paralogs and those that are not paralogs...')
f.write('%g, %g \n' %(stats.mannwhitneyu(ncons.pident[ncons.qseqid.isin(par)], ncons1.pident)))
f.write('Two sample wilcoxon test between percentage identity of conserved genes that are paralogs and those that are not paralogs...')
f.write('%g, %g \n' %(stats.mannwhitneyu(cons.pident[cons.qseqid.isin(par)], cons1.pident)))
f.write('Two sample wilcoxon test between cgo non paralogs and ncgo non paralogs')
f.write('%g, %g \n' %(stats.mannwhitneyu(ncons1.pident, cons1.pident)))
f.close()
#cons.qseqid[cons.qseqid.isin(par)].to_csv('cons_par1.txt', sep='\t',index=None)
#ncons.qseqid[ncons.qseqid.isin(par)].to_csv('ncons_par1.txt', sep='\t',index=None)


sns.palplot(sns.color_palette("magma",n_colors=4))

a4_dims = (10, 8)
sns.set_style('white')
fig, ax = pyplot.subplots(figsize=a4_dims)
sns.set_palette(palette="magma",n_colors=4)
sns.distplot(cons.pident[~cons.qseqid.isin(par)], hist=False, rug=True, kde_kws=dict(linewidth=1.25), label="CGOs in single copy")
sns.distplot(cons.pident[cons.qseqid.isin(par)], hist=False, rug=True,kde_kws=dict(linewidth=1.25),label="CGOs with paralog")
sns.distplot(ncons.pident[~ncons.qseqid.isin(par)], hist=False, rug=True,kde_kws=dict(linewidth=1.25),label="non-CGOs in single copy")
sns.distplot(ncons.pident[ncons.qseqid.isin(par)], hist=False, rug=True,kde_kws=dict(linewidth=1.25),label="non-CGOs with paralog")

ax.set_title('Distriburion plot')
ax.set_xlabel('Percentage Identity')

# m-m group
sns.palplot(sns.color_palette("magma",n_colors=4))
a4_dims = (10, 8)
sns.set_style('white')
fig, ax = pyplot.subplots(figsize=a4_dims)
sns.set_palette(palette="magma",n_colors=2)
#sns.distplot(cons.pident, hist=False, rug=True, kde_kws=dict(linewidth=1.25), label="CGOs")
#sns.distplot(cons.pident[cons.qseqid.isin(par)], hist=False, rug=True,kde_kws=dict(linewidth=1.25),label="CGOs with Paralog")
sns.distplot(identity_cgo, hist=False, rug=True,kde_kws=dict(linewidth=1.25),label="CGOs")
sns.distplot(identity_ncgo, hist=False, rug=True,kde_kws=dict(linewidth=1.25),label="Non-CGOs")
ax.set_title('Distriburion plot')
ax.set_xlabel('Percentage Identity')

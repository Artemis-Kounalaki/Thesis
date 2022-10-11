import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os


pd.set_option('display.max_rows', 50)
pd.set_option('display.max_columns', 50)

# No hits in macaca human - macaca
os.chdir(os.path.expanduser('~/conserved_gene_order1/macaca_reference'))
df_h_mac = pd.read_csv('CGO_nCGO_h-mac.txt', sep='\t')
hits_mac_h = df_h_mac.iloc[:, 1].tolist()


ids_macaca = pd.read_csv('clean_ids_macaca.txt', sep='\t')
ids_macaca = ids_macaca.iloc[:, 1]


no_hits_mac_h = ids_macaca[~ids_macaca.isin(hits_mac_h)]


# No hits in mus human-mus
os.chdir(os.path.expanduser('~/conserved_gene_order1/mus_reference'))
df_h_mus = pd.read_csv('CGO_nCGO_h-mus.txt', sep='\t')
hits_mus_h = df_h_mus.iloc[:, 1].tolist()


ids_mus = pd.read_csv('clean_ids_mus.txt', sep='\t')
ids_mus = ids_mus.iloc[:, 1]


no_hits_mus_h = ids_mus[~ids_mus.isin(hits_mus_h)]


# No hits mus-macaca
# mac
os.chdir(os.path.expanduser('~/conserved_gene_order1/mus_reference'))
df_m = pd.read_csv('CGO_nCGO_m-m.txt', sep='\t')
hits_m_mac = df_m.iloc[:, 1].tolist()
no_hits_mac_m = ids_macaca[~ids_macaca.isin(hits_m_mac)]

# mus
hits_m_mus = df_m.iloc[:, 0].tolist()
no_hits_mus_m = ids_mus[~ids_mus.isin(hits_m_mus)]


# No hits common in all groups

# No hits in human human - macaca
os.chdir(os.path.expanduser('~/conserved_gene_order1/human_reference'))
hits_hum_mac = df_h_mac.iloc[:, 0].tolist()

ids_hum = pd.read_csv('clean_ids.txt', sep='\t')
ids_hum = ids_hum.iloc[:, 1]


no_hits_hum_mac = ids_hum[~ids_hum.isin(hits_hum_mac)]


# No hits in human human - mus
hits_hum_mus = df_h_mus.iloc[:, 0].tolist()


no_hits_hum_mus = ids_hum[~ids_hum.isin(hits_hum_mus)]

common = pd.merge(no_hits_hum_mus, no_hits_hum_mac, how='inner')

# Human Paralogs in no hits


# Find paralogs

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
df_paralogs["pairs"] = df_paralogs.apply(lambda row: ''.join(sorted([row["qseqid"], row["sseqid"]])), axis=1)
only_pairs = df_paralogs[df_paralogs["pairs"].duplicated(keep = False)].sort_values(by = "pairs") # find reciprocallity
keep=list(set(only_pairs['qseqid']).intersection(set(only_pairs['sseqid']))) #prevent false pairs ,alphabetically sort pairs
only_pairs=only_pairs[only_pairs['qseqid'].isin(keep) & only_pairs['sseqid'].isin(keep)]
only_pairs = only_pairs.reset_index(drop=True)
par = pd.unique(only_pairs[['qseqid', 'sseqid']].values.ravel())


# Plot no hits in barplot

x = [u'Macaca Proteins in H-Mac', u'Macaca  proteins in Mus-Macaca', u'Common no hits',u'Mus proteins in H-MUs',u'Mus proteins in Mus-Mac',u'Common no hits', u'Human proteins in H-Mus',u'Human proteins in H-Macaca',u'Common no hits',u'Common no hits with Paralog']
y = [len(no_hits_mac_h), len(no_hits_mac_m), len(pd.merge(no_hits_mac_h, no_hits_mac_m, how='inner')),len(no_hits_mus_h), len(no_hits_mus_m),len(pd.merge(no_hits_mus_h, no_hits_mus_m, how='inner')),len(no_hits_hum_mus),len(no_hits_hum_mac),len(pd.merge(no_hits_hum_mus, no_hits_hum_mac, how='inner')),len(common[common.ID.isin(par)])]
fig, ax = plt.subplots()
width = 0.5
ind = np.arange(len(y))
ax.barh(ind, y, width, color="navajowhite")
ax.set_yticks(ind+width/2)
ax.set_yticklabels(x, minor=False)
plt.title('No hits')
plt.xlabel('Number of Genes')
plt.ylabel('Groups')
for i, v in enumerate(y):
    ax.text(v +5, i - .1, str(v), color='indianred', fontweight='bold')
plt.savefig(os.path.join('nohits_plot.png'), dpi=300, format='png', bbox_inches='tight')


# 1391 found as common paralogs analysis

com = common[common.ID.isin(par)]
com1 = pd.unique(com.values.ravel())
matches = only_pairs[only_pairs['qseqid'].isin(com1)]
matches2 = only_pairs[only_pairs['sseqid'].isin(com1)]
s1 = pd.merge(matches, matches2, how='inner', on=['qseqid', 'sseqid'])
s1.rename(columns={'pairs_x': 'pairs'}, inplace=True)

# Create a df with common no hits and their paralogs

spec1 = only_pairs[['qseqid', 'sseqid']][only_pairs['qseqid'].isin(com1)]
sepc1 = spec1.sort_values(['qseqid'], ascending=[True])
spec1['CNH'] = (spec1.sseqid.isin(com1)).astype('int')
spec1=spec1.groupby(['qseqid']).agg({'sseqid': list, 'CNH': sum}).reset_index()
spec1['matches']=spec1.sseqid.str.len()


# Same matches as no hits

g1=spec1[spec1.matches==spec1.CNH]
g2=spec1[spec1.matches-spec1.CNH==1]
g3=spec1[spec1.matches-spec1.CNH>1]
g4=spec1[spec1.CNH==0]
g1.to_csv('~/CNHall.txt', header=None, index=None)
g2.to_csv('~/CNH_1hit.txt', header=None, index=None)
g3.to_csv('~/CNH>1hit.txt', header=None, index=None)
g4.to_csv('~/1CNH.txt', header=None, index=None)

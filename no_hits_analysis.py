import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# No hits in macaca human - macaca
df_h_mac = pd.read_csv('al.txt', sep='\t')
hits_mac_h= df_h_mac.iloc[:,1].tolist()


ids_macaca= pd.read_csv('clean_ids_macaca.txt', sep='\t')
ids_macaca = ids_macaca.iloc[:,1]


no_hits_mac_h= ids_macaca[~ids_macaca.isin(hits_mac_h)]


# No hits in mus human-mus

df_h_mus = pd.read_csv('al_mus.txt', sep='\t')
hits_mus_h= df_h_mus.iloc[:,1].tolist()


ids_mus= pd.read_csv('clean_ids_mus.txt', sep='\t')
ids_mus = ids_mus.iloc[:,1]


no_hits_mus_h= ids_mus[~ids_mus.isin(hits_mus_h)]



# No hits mus-macaca
# mac
df_m = pd.read_csv('al_m-m.txt', sep='\t')
hits_m_mac= df_m.iloc[:,1].tolist()
no_hits_mac_m = ids_macaca[~ids_macaca.isin(hits_m_mac)]

# mus
hits_m_mus= df_m.iloc[:,0].tolist()
no_hits_mus_m = ids_mus[~ids_mus.isin(hits_m_mus)]


# No hits common in all groups

print(pd.merge(no_hits_mac_h, no_hits_mac_m, how='inner'))
print(pd.merge(no_hits_mus_h, no_hits_mus_m, how='inner'))


# No hits in human human - macaca

df_h_mac = pd.read_csv('al.txt', sep='\t')
hits_hum_mac= df_h_mac.iloc[:,0].tolist()

ids_hum= pd.read_csv('clean_ids_human.txt', sep='\t')
ids_hum = ids_hum.iloc[:,1]
print(ids_hum)

no_hits_hum_mac= ids_hum[~ids_hum.isin(hits_hum_mac)]
print(no_hits_hum_mac)

# No hits in human human - mus
df_h_mus = pd.read_csv('al_mus.txt', sep='\t')
hits_hum_mus= df_h_mus.iloc[:,0].tolist()


no_hits_hum_mus= ids_hum[~ids_hum.isin(hits_hum_mus)]
print(no_hits_hum_mus)

common=pd.merge(no_hits_hum_mus, no_hits_hum_mac, how='inner')

# HUman Paralogs in no hits


# Find paralogs

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
par= pd.unique(only_pairs[['qseqid','sseqid']].values.ravel())

print(common[common.ID.isin(par)])


x = [u'Macaca in H-Mac', u'Macaca in Mus-Macaca', u'Common',u'Mus in H-MUs',u'Mus in Mus-Mac',u'Common', u'Human in H-Mus',u'Human in H-Macaca',u'Common',u'Common Paralogs']
y = [len(no_hits_mac_h), len(no_hits_mac_m), len(pd.merge(no_hits_mac_h, no_hits_mac_m, how='inner')),len(no_hits_mus_h), len(no_hits_mus_m),len(pd.merge(no_hits_mus_h, no_hits_mus_m, how='inner')),len(no_hits_hum_mus),len(no_hits_hum_mac),len(pd.merge(no_hits_hum_mus, no_hits_hum_mac, how='inner')),len(common[common.ID.isin(par)])]
fig, ax = plt.subplots()
width = 0.5
ind = np.arange(len(y))
ax.barh(ind, y, width, color="navajowhite")
ax.set_yticks(ind+width/2)
ax.set_yticklabels(x, minor=False)
plt.title('No hits')
plt.xlabel('Genes')
plt.ylabel('Groups')
for i, v in enumerate(y):
    ax.text(v +5, i - .1, str(v), color='red', fontweight='bold')
plt.savefig(os.path.join('test.png'), dpi=300, format='png', bbox_inches='tight')

# There are four groups derived from hits in its group.
# In its group how many genes from each chromosome belong to these subgroups?

import pandas as pd
from scipy import stats
from scipy.stats import mannwhitneyu
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


# Find paralogs
print('Loading human paralogs ...')

#os.chdir(os.path.expanduser(path_results_paralogs))
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


# Statistic paralog conserved and paralog non conserved

# Loading results al.txt
df_all = pd.read_csv('CGO_nCGO_h-mac.txt', sep='\t')
df_all.columns=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'slen', 'qlen', 'Conserved']
genes=df_all[['qseqid','sseqid']].values.ravel()
par= pd.unique(only_pairs[['qseqid','sseqid']].values.ravel())

df_all['Paralog']= df_all.qseqid.isin(par)
df_all['Paralog']= df_all['Paralog'].fillna(0)
df_all['Paralog']= df_all['Paralog'].astype(int)
cons= df_all[df_all.Conserved==1]
cons.pident = cons.pident.astype(float)
ncons= df_all[df_all.Conserved==0]
ncons.pident = ncons.pident.astype(float)
cons1=cons[~cons.qseqid.isin(par)]
ncons1=ncons[~ncons.qseqid.isin(par)]

# 4 groups
# GCOs & Paralogous

print('The number of conserved genes that are paralogs is:')
print(len(cons[cons.qseqid.isin(par)]))
print('Their average percentage identity is:',sum(cons.pident[cons.qseqid.isin(par)])/len(cons.pident[cons.qseqid.isin(par)]))

# nGCOs & Paralogous

print('The number of non conserved genes that are paralogs is:')
print(len(ncons[ncons.qseqid.isin(par)]))
print('Their average percentage isentity is:',sum(ncons.pident[ncons.qseqid.isin(par)])/len(ncons.pident[ncons.qseqid.isin(par)]))

# GCOs & Non Paralogous

print('The number of conserved genes without paralogs is:',len(cons1))
print('The average percentage identity of conserved genes without paralogs is:',sum(cons1.pident)/len(cons1))

# nGCOs & Non Paralogous

print('The number of non conserved genes without paralogs is:',len(ncons1))


# Loading human ids

df_human = pd.read_csv('clean_ids.txt', sep='\t', index_col=0)
hits=pd.unique(df_all.qseqid)
df_nohits = df_human[~df_human.ID.isin(hits)]

CGO=df_human[df_human.ID.isin(cons.qseqid)]
nCGO=df_human[df_human.ID.isin(ncons.qseqid)]
CGO_Par = df_human[df_human.ID.isin(cons.qseqid[cons.qseqid.isin(par)])]
nCGO_Par = df_human[df_human.ID.isin(ncons.qseqid[ncons.qseqid.isin(par)])]
nCGO_nPar= df_human[df_human.ID.isin(ncons.qseqid[~ncons.qseqid.isin(par)])]
CGO_nPar = df_human[df_human.ID.isin(cons.qseqid[~cons.qseqid.isin(par)])]
CGO_Par = CGO_Par.reset_index(drop=True)
nCGO_Par = nCGO_Par.reset_index(drop=True)
nCGO_nPar = nCGO_nPar.reset_index(drop=True)
CGO_nPar = CGO_nPar.reset_index(drop=True)

# How many GCOs paralogs have GCOs & nGCOs

df_par=pd.DataFrame()
df_par['qseqid']=par[::2]
df_par['sseqid']=par[1::2]


one=df_par.index[df_par['qseqid'].isin(CGO_Par.ID)].tolist()
two=df_par.index[df_par['sseqid'].isin(nCGO_Par.ID)].tolist()
c = set(one) & set(two)


df_c_n=pd.DataFrame(columns=['CGO','nCGO'])

count=0
for i in c:
    count+=1
    if (df_par['qseqid'].iloc[i] in CGO_Par.ID.tolist()):
        df_c_n.loc[count,'CGO']= df_par['qseqid'].iloc[i]
        df_c_n.loc[count,'nCGO']= df_par['sseqid'].iloc[i]
    else:
        df_c_n.loc[count,'CGO']= df_par['sseqid'].iloc[i]
        df_c_n.loc[count,'nCGO']= df_par['qseqid'].iloc[i]

#df_c_n.to_csv('/Users/artemiskounalake/cgo_ncgo_par.txt', header=None, index=None, sep=' ', mode='a')

# end here


# Code for Chromosome plot

# Chromosomes

list_chrom=range(1,23)
len_chr_g=[]
len_chr=[]


for i in list_chrom:
    len_chr_g.append(len(df_human[df_human.Chrom==str(i)]))
len_chr_g.append(len(df_human[df_human.Chrom=='X']))
len_chr_g.append(len(df_human[df_human.Chrom=='Y']))


def Chrom(lis):
    lis_c=[]
    for i in list_chrom:
        lis_c.append(len(lis[lis['Chrom']==str(i)]))
    lis_c.append(len(lis[lis['Chrom']=='X']))
    lis_c.append(len(lis[lis['Chrom']=='Y']))
    return lis_c

#One chromosome plot of subgroups

CGO_Par = np.array(np.array(Chrom(CGO_Par))/np.array(len_chr_g))*100
CGO_nPar = np.array(np.array(Chrom(CGO_nPar))/np.array(len_chr_g))*100
nCGO_Par = np.array(np.array(Chrom(nCGO_Par))/np.array(len_chr_g))*100
nCGO_nPar = np.array(np.array(Chrom(nCGO_nPar))/np.array(len_chr_g))*100
df_nohits = np.array(np.array(Chrom(df_nohits))/np.array(len_chr_g))*100
x=[u'1', u'2', u'3', u'4', u'5', u'6', u'7', u'8', u'9', u'10', u'11', u'12', u'13', u'14', u'15', u'16', u'17', u'18', u'19', u'20', u'21', u'22',  u'X', u'Y']
chroms=['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22',  'X', 'Y']
fig = plt.figure(figsize=(11,10), dpi=200)
left, bottom, width, height = 0.2, 0.1, 0.8, 0.5
ax = fig.add_axes([left, bottom, width, height])
width = 0.6
ax.set_xlabel('Percentage in groups')
ax.set_yticks(np.arange(len(chroms)))
ax.set_yticklabels(chroms)
ax.set_ylabel('Chromosome')

ticks = range(len(chroms))
ax.barh(ticks, CGO_Par, width, left=CGO_nPar+nCGO_Par+nCGO_nPar+df_nohits, color='thistle',label='CGOs with Par')
ax.barh(ticks, CGO_nPar, width, align='center', left=nCGO_Par+nCGO_nPar+df_nohits,color='powderblue',label='CGOs without Par')
ax.barh(ticks, nCGO_Par, width, align="center",left=nCGO_nPar+df_nohits,color='bisque',label='nCGOs with Par')
ax.barh(ticks, nCGO_nPar, width, align="center",left=df_nohits,color='orchid',label='nCGOs without Par')
ax.barh(ticks, df_nohits, width, align="center",color='silver',label='No Hits')
ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.10),ncol=3, fancybox=True, shadow=True)

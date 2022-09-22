import pandas as pd
import numpy as np
from scipy import stats

# Find best reciprocal couple of paralogs in human.

df_paralogs = pd.read_csv('clean_blast_results_human.txt', sep=' ')
print(df_paralogs)
'''
df_paralogs.columns=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'slen', 'qlen']
df_paralogs.qseqid = df_paralogs.qseqid .astype(str)
df_paralogs.sseqid = df_paralogs.sseqid .astype(str)
df_del = df_paralogs[df_paralogs.qseqid == df_paralogs.sseqid] # delete the hit of a protein with itself
df_paralogs=df_paralogs[~df_paralogs.isin(df_del)] 
df_paralogs = df_paralogs.dropna(how='all')  
df_paralogs = df_paralogs.reset_index(drop=True) 
df_paralogs["pairs"] = df_paralogs.apply(lambda row: ''.join(sorted([row["qseqid"], row["sseqid"]])), axis = 1)
only_pairs = df_paralogs[df_paralogs["pairs"].duplicated(keep = False)].sort_values(by = "pairs") # find reciprocallity
only_pairs=only_pairs[~only_pairs.qseqid.duplicated(keep=False)] 

only_pairs = only_pairs.reset_index(drop=True) 

print(only_pairs)
# Now upload the results of best reciprocal blast (Human-Macaca) with information : Conserved/non Conserved

df_all = pd.read_csv('al0.txt', sep='\t')
df_all.columns=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'slen', 'qlen', 'Conserved']



# Place to results the paralog results of human


df_all['Paralog'] = (df_all.qseqid).map(only_pairs.set_index('qseqid')['sseqid']) # add paralog of human protein
df1= df_all[['Paralog','sseqid']] # keep the couples (paralog-macaca protein)
df1.columns=['qseqid', 'sseqid']
df_all['P_ident']= (df1.qseqid).map(df_all.set_index('qseqid')['pident'])
df_all['Cons']= (df1.qseqid).map(df_all.set_index('qseqid')['Conserved'])
'''
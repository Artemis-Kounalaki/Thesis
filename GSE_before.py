import pandas as pd
from scipy import stats

# Groups preparation for GSE study

# Load RBH results with GCO info (1:GCO,0:NGCO)
df_al = pd.read_csv('CGO_nCGO_h-mus.txt', sep='\t')
df_al.columns=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'slen', 'qlen','Conserved']
print(df_al)

# Keep only human proteins have same hit on macaca
df_dou = df_al[df_al.duplicated(subset=['sseqid'], keep=False)]
df_dou= df_dou.reset_index(drop=True)
print(df_dou.sort_values('sseqid'))


# Cons/ non Cons

one=df_dou[df_dou.Conserved == 1]
two=df_dou[df_dou.Conserved == 0]
print(sum(one.pident)/len(one))
print(sum(two.pident)/len(two))
print(stats.mannwhitneyu(one.pident, two.pident))

# Their average pr identity & pvalue
target_c = df_dou[df_dou.Conserved == 1]
target_nc = df_dou[df_dou.Conserved == 0]

# Total num of cons
print(target_c)

# Total num of NON cons
print(target_nc)


# Non uniform group number of groups
keep=set(one['sseqid']).intersection(two['sseqid'])
print(len(keep))

# Total number of non uniform prot
inte=df_dou[df_dou.sseqid.isin(keep)]

# Num of con in non uni
c_nuni=inte[inte.Conserved ==1].qseqid
print(c_nuni)

# Num of non cons in non uni
nc_nuni=inte[inte.Conserved ==0].qseqid
print(nc_nuni)


# Total number  in uniform prot
res=df_dou[~df_dou.sseqid.isin(keep)]
print(res)


# Groups of cons in uniform
print(len(set(res[res.Conserved ==1].sseqid)))
# Groups of non cons in uniform
print(len(set(res[res.Conserved ==0].sseqid)))


# Num of con in uni
c_uni=res[res.Conserved ==1].qseqid
print(c_uni)

# Num of non cons in uni
nc_uni=res[res.Conserved ==0].qseqid
print(nc_uni)

# Save cons/ncons non uni group
c_nuni.to_csv('~/c_nuni.txt', header=None, index=None)
nc_nuni.to_csv('~/nc_nuni.txt', header=None, index=None)

#Save cons/ncons uni groups
c_uni.to_csv('~/c_uni.txt', header=None, index=None)
nc_uni.to_csv('~/nc_uni.txt', header=None, index=None)

#Save cons/ncons
target_c.qseqid.to_csv('~/target_c.txt', header=None, index=None)
target_nc.qseqid.to_csv('~/target_nc.txt', header=None, index=None)

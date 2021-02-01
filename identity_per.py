import pandas as pd
import numpy as np

os.chdir(os.path.expanduser('~/conserved_gene_order/macaca_reference'))
cgo = pd.read_csv('CGO.txt', sep=" ", header=None)
cgo.columns = ["Prot_human", "Prot_macaca"]
cgo = cgo.sort_values(by=['Prot_human'])
cgo = cgo.reset_index(drop=True)

ncgo = pd.read_csv('nCGO.txt', sep=" ", header=None)
ncgo.columns = ["Prot_human", "Prot_macaca"]
ncgo = ncgo.sort_values(by=['Prot_human'])
ncgo = ncgo.reset_index(drop=True)

reci = pd.read_csv('Reciprocal_h-m.txt', sep=" ", header=None)
reci.columns = ["Prot_human", "Prot_macaca","P.identity", 'Length', '4', 'Gap', '6', '7', '8', '9', '10', '11', 'Hum_len', 'Mac_len']


Cgo = pd.merge(reci, cgo, on=['Prot_human','Prot_macaca'], how='left', indicator='Exist')
Cgo['Exist'] = np.where(Cgo.Exist == 'both', True, False)

Cgo=Cgo.loc[Cgo['Exist'] == True]
Cgo=Cgo.drop(columns=['4', '6', '7', '8', '9', '10', '11'])
Cgo = Cgo.drop_duplicates(subset = ["Prot_human",'Prot_macaca'])
Cgo = Cgo.reset_index(drop=True)
Cgo['Conserved']='1'
nCgo = pd.merge(reci, ncgo, on=['Prot_human','Prot_macaca'], how='left', indicator='Exist')
nCgo['Exist'] = np.where(nCgo.Exist == 'both', True, False)

nCgo=nCgo.loc[nCgo['Exist'] == True]
nCgo=nCgo.drop(columns=[ '4','6', '7', '8', '9', '10', '11'])
nCgo = nCgo.drop_duplicates(subset = ["Prot_human",'Prot_macaca'])
nCgo = nCgo.reset_index(drop=True)
nCgo['Conserved']='0'

df_al=pd.concat([Cgo, nCgo], ignore_index=True)
df_al=df_al.drop(columns=['Exist'])
df_al = df_al.reset_index(drop=True)


print(Cgo['P.identity'].sum()/len(Cgo))
print(nCgo['P.identity'].sum()/len(nCgo))

df_al.to_csv('al.txt', sep='\t')

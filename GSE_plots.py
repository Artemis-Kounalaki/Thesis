# libraries & dataset
import pandas as pd
import seaborn as sns
import numpy as np
from scipy.stats import zscore
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from sklearn.preprocessing import LabelEncoder

# Load r results

df_expr = pd.read_csv('data_expr.txt', sep=' ')
cgo=pd.read_csv('c_nuni_2971.txt', sep=',', header=None)
cgo.columns=['name']
cgo = cgo.dropna(axis=0)
ncgo=pd.read_csv('nc_nuni_2971.txt', sep=',', header=None)
ncgo.columns=['name']


df_expr = df_expr[df_expr.index.notna()]
df_expr.reset_index(level=0, inplace=True)
df_expr['name']= df_expr['index']
df_expr=df_expr.drop(columns=['index'])

out1 = (cgo.merge(df_expr, right_on='name',left_on='name'))
print(out1)
out1=out1.drop_duplicates(subset='name', keep="last")
out1.replace([np.inf, -np.inf], np.nan, inplace=True)
print(out1)

out2 = (ncgo.merge(df_expr, right_on='name',left_on='name'))
print(out2)
out2=out2.drop_duplicates(subset='name', keep="last")
out2.replace([np.inf, -np.inf], np.nan, inplace=True)
print(out2)

#for GC step
#out1.name.to_csv('~/target_c.txt', header=None, index=None)
#out2.name.to_csv('~/target_nc.txt', header=None, index=None)

# HISTPLOTS in 1
names = list(out1.columns.values)
names.remove('name')
suffix = 'Normal human tissue:'
names2 = [x[len(suffix):] for x in names if x.startswith(suffix)]
count=-1
fig, axs = plt.subplots(6,6, figsize=(20, 20))
for i in range(0,6):
    for y in range(0,6):
        count+=1
        sns.histplot(out1[names[count]], kde=True,color='plum',ax=axs[i,y]).set_label(" ")
        sns.histplot(out2[names[count]], kde=True,color='skyblue', ax=axs[i,y]).set_label(" ")
        plt.legend(["GCOs","nGCOs"],bbox_to_anchor = (1.05, 0.6))


#CLUSTERMAP/PCA
out1.insert(0, 'CGO', '1')
out2.insert(0, 'CGO', '0')
out1=out1.T
out2=out2.T
out1.columns=out1.iloc[1,:].tolist()
out1.drop("name",axis=0,inplace=True)
out2.columns=out2.iloc[1,:].tolist()
out2.drop("name",axis=0,inplace=True)
tab = pd.concat([out1, out2], axis=1, join="inner")
tab=tab.T
tab=tab.fillna(0)
tab = tab.sort_values(by ='CGO' )
target=tab['CGO']
my_palette = {'0' : 'lightcoral', '1' :'blue'}
row_colors = target.map(my_palette)
tab=tab.drop("CGO",axis=1)
tab=tab.apply(zscore)
sns.set(font_scale=0.75)

#CLUSTERMAP
sns.clustermap(tab,metric='correlation', method="average", row_colors = row_colors, row_cluster = False, cmap='pink',figsize=(15,10))
#sns_plot1.savefig("heatmap01.png")
sns.clustermap(tab,metric='correlation', method="average", row_colors = row_colors, col_cluster = False, cmap='pink',figsize=(12,18))
#sns_plot2.savefig("heatmap10.png")
sns.clustermap(tab,metric='correlation', method="average", row_colors = row_colors, standard_scale=0, cmap='pink',figsize=(10,10))
#sns_plot3.savefig("heatmap11-0.png")
sns.clustermap(tab,metric='correlation', method="average", row_colors = row_colors, standard_scale=1, cmap='pink',figsize=(10,10))
#sns_plot4.savefig("heatmap11-1.png")
plt.show()



# PCA

pca = PCA(n_components=36)
X = pca.fit(tab).transform(tab)
le = LabelEncoder()
le.fit(target)
x_lan=le.transform(target)
pca_df = pd.DataFrame(columns = ['x', 'y', 'name', 'label'])
print(pca_df)
pca_df['PCA1'] = X[:, 0]
pca_df['PCA2'] = X[:, 1]
pca_df['label'] = x_lan
sns.set(style='whitegrid', palette='PuRd')
ax = sns.scatterplot(x='PCA1', y='PCA2', hue='label', data=pca_df)
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.show()

import pandas as pd
import seaborn as sns
from matplotlib import pyplot
from scipy.stats import zscore
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt

# Load r results
df_expr = pd.read_csv('data_expr.txt', sep=' ')
cgo=pd.read_csv('target_cgo.txt', sep=',', header=None)
cgo.columns=['name']
ncgo=pd.read_csv('target_ncgo.txt', sep=',', header=None)
ncgo.columns=['name']

df_expr = df_expr[df_expr.index.notna()]
df_expr.reset_index(level=0, inplace=True)
df_expr['name']= df_expr['index']
df_expr=df_expr.drop(columns=['index'])

out1 = (cgo.merge(df_expr, right_on='name',left_on='name'))
out1=out1.drop_duplicates(subset='name', keep="last")


out2 = (ncgo.merge(df_expr, right_on='name',left_on='name'))
out2=out2.drop_duplicates(subset='name', keep="last")



sns.palplot(sns.color_palette("ch:s=-.2,r=.6"))
names = list(out1.columns.values)
names.remove('name')

for name in names:
    a4_dims = (5, 2)
    sns.set_style('white')
    fig, ax = pyplot.subplots(figsize=a4_dims)
    sns.set_palette("ch:s=-.2,r=.6")
    #sns.distplot(cons.pident, hist=False, rug=True, kde_kws=dict(linewidth=1.25), label="CGOs")
    #sns.distplot(cons.pident[cons.qseqid.isin(par)], hist=False, rug=True,kde_kws=dict(linewidth=1.25),label="CGOs with Paralog")
    sns.distplot(out1[name], hist=False, rug=True,kde_kws=dict(linewidth=1.25),label="CGOs with best paralog in non-CGOs")
    sns.distplot(out2[name], hist=False, rug=True,kde_kws=dict(linewidth=1.25),label="non-CGOs with best paralog in CGOs")
    ax.set_title('Distriburion plot in %s'%name)
    ax.set_xlabel('Gene expression')

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

# Prepare a vector of color mapped to the 'cyl' column
target=tab['CGO']
my_palette = {'0' : 'lightcoral', '1' :'pink'}
row_colors = target.map(my_palette)
print(row_colors)
tab=tab.drop("CGO",axis=1)
tab=tab.apply(zscore)

#tab['CGO'] = tab['CGO'].astype(float)
sns.set(font_scale=0.75)

#sns.clustermap(tab, metric="correlation", method="single", cmap="PuRd",  row_colors=row_colors)
sns_plot1=sns.clustermap(tab,metric='correlation', method="average", row_colors = row_colors, row_cluster = False, cmap='pink',figsize=(15,10))
#sns_plot1.savefig("heatmap01.png")
sns_plot2=sns.clustermap(tab,metric='correlation', method="average", row_colors = row_colors, col_cluster = False, cmap='pink',figsize=(12,18))
#sns_plot2.savefig("heatmap10.png")
sns_plot3=sns.clustermap(tab,metric='correlation', method="average", row_colors = row_colors, standard_scale=0, cmap='pink',figsize=(10,10))
#sns_plot3.savefig("heatmap11-0.png")
sns_plot4=sns.clustermap(tab,metric='correlation', method="average", row_colors = row_colors, standard_scale=1, cmap='pink',figsize=(10,10))
#sns_plot4.savefig("heatmap11-1.png")

plt.show()

tab=tab.drop(tab.var()[(tab.var() <0.3)].index, axis=1)


# PCA


pca = PCA(n_components=36)
X = pca.fit(tab).transform(tab)
from sklearn.preprocessing import LabelEncoder
le = LabelEncoder()
le.fit(target)
x_lan=le.transform(target)

pca_df = pd.DataFrame(columns = ['x', 'y', 'name', 'label'])
pca_df['PCA1'] = X[:, 0]
pca_df['PCA2'] = X[:, 1]
pca_df['CGO'] = target
pca_df['label'] = x_lan
sns.set(style='whitegrid', palette='magma')
ax = sns.scatterplot(x='PCA1', y='PCA2', hue='label', data=pca_df)
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.show()


# Kmeans

pca = PCA(n_components=2)
tab = StandardScaler().fit_transform(tab)
principalComponents = pca.fit_transform(tab)
kmeans = KMeans(n_clusters=2, max_iter=300, algorithm = 'auto').fit(principalComponents)
labels=kmeans.labels_

for i in range(0, principalComponents.shape[0]):
    if labels[i] == 0:
        c1 = plt.scatter(principalComponents[i,0],principalComponents[i,1],c='violet',marker='+')
    elif labels[i] == 1:
        c2 = plt.scatter(principalComponents[i,0],principalComponents[i,1],c='c',marker='+')

plt.scatter(kmeans.cluster_centers_[:, 0], kmeans.cluster_centers_[:, 1],marker = "d", c='b', s=50)
plt.legend()
plt.title('K means clusters dataset into' + str(2) +  'clusters')
plt.show()

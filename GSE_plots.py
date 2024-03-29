import pandas as pd
import seaborn as sns
import numpy as np
from scipy.stats import zscore
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from sklearn.preprocessing import LabelEncoder


# Load r results

df_expr = pd.read_csv('data_expr_2971.txt', sep=' ')
print(df_expr)
c_uni=pd.read_csv('c_uni_2971.txt', sep=',', header=None)
c_uni.columns=['name']
c_uni = c_uni.dropna(axis=0)

nc_uni=pd.read_csv('nc_uni_2971.txt', sep=',', header=None)
nc_uni.columns=['name']

c_nuni=pd.read_csv('c_nuni_2971.txt', sep=',', header=None)
c_nuni.columns=['name']

nc_nuni=pd.read_csv('nc_nuni_2971.txt', sep=',', header=None)
nc_nuni.columns=['name']


df_expr = df_expr[df_expr.index.notna()]
df_expr.reset_index(level=0, inplace=True)
df_expr['name']= df_expr['index']
df_expr=df_expr.drop(columns=['index'])
print(df_expr)

out1 = (c_uni.merge(df_expr, right_on='name',left_on='name'))
out1=out1.drop_duplicates(subset='name', keep="last")
out1.replace([np.inf, -np.inf], np.nan, inplace=True)
print(out1)

out2 = (nc_uni.merge(df_expr, right_on='name',left_on='name'))
out2=out2.drop_duplicates(subset='name', keep="last")
out2.replace([np.inf, -np.inf], np.nan, inplace=True)
print(out2)

out3 = (c_nuni.merge(df_expr, right_on='name',left_on='name'))
out3=out3.drop_duplicates(subset='name', keep="last")
out3.replace([np.inf, -np.inf], np.nan, inplace=True)
print(out3)

out4 = (nc_nuni.merge(df_expr, right_on='name',left_on='name'))
out4=out4.drop_duplicates(subset='name', keep="last")
out4.replace([np.inf, -np.inf], np.nan, inplace=True)
print(out4)


# for GC step
#out1.name.to_csv('~/target_c_3416.txt', header=None, index=None)
#out2.name.to_csv('~/target_nc_3416.txt', header=None, index=None)


names = list(out1.columns.values)
names.remove('name')
suffix = 'Normal human tissue:'
names2 = [x[len(suffix):] for x in names if x.startswith(suffix)]

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
        sns.histplot(out1[names[count]], kde=True,color='red',ax=axs[i,y]).set_label(" ")
        sns.histplot(out2[names[count]], kde=True,color='skyblue', ax=axs[i,y]).set_label(" ")
        sns.histplot(out3[names[count]],kde=True,color='gray',ax=axs[i,y]).set_label(" ")
        sns.histplot(out4[names[count]],kde=True,color='plum', ax=axs[i,y]).set_label(" ")
        plt.legend(["GCOs in Uniform","nGCOs in Uniform","GCOs in non-Uniform ","nGCOs in non-Uniform"],bbox_to_anchor = (1.05, 0.6))

out1.insert(0, 'CGO', '1')
out2.insert(0, 'CGO', '0')
out3.insert(0, 'CGO', '2')
out4.insert(0, 'CGO', '3')


out1=out1.T
out2=out2.T
out3=out3.T
out4=out4.T

out1.columns=out1.iloc[1,:].tolist()
out1.drop("name",axis=0,inplace=True)

out2.columns=out2.iloc[1,:].tolist()
out2.drop("name",axis=0,inplace=True)

out3.columns=out3.iloc[1,:].tolist()
out3.drop("name",axis=0,inplace=True)

out4.columns=out4.iloc[1,:].tolist()
out4.drop("name",axis=0,inplace=True)
print(out1,out2,out3,out4)


frames=[out1,out2,out3,out4]
tab = pd.concat(frames, axis=1, join="inner")
tab=tab.T
tab=tab.fillna(0)
tab = tab.sort_values(by ='CGO' )


# Prepare a vector of color mapped to the 'cyl' column
target=tab['CGO']
my_palette = {'0' : 'skyblue', '1' :'red','2' : 'gray', '3' :'plum'}
row_colors = target.map(my_palette)
tab=tab.drop("CGO",axis=1)
tab=tab.apply(zscore)
sns.set(font_scale=0.75)
sns_plot1=sns.clustermap(tab,metric='correlation', method="average", row_colors = row_colors, row_cluster = False, cmap='pink',figsize=(15,10))
sns_plot1.savefig("heatmap01.png")
sns_plot2=sns.clustermap(tab,metric='correlation', method="average", row_colors = row_colors, col_cluster = False, cmap='pink',figsize=(12,18))
sns_plot2.savefig("heatmap10.png")
sns_plot3=sns.clustermap(tab,metric='correlation', method="average", row_colors = row_colors, standard_scale=0, cmap='pink',figsize=(10,10))
sns_plot3.savefig("heatmap11-0.png")
sns_plot4=sns.clustermap(tab,metric='correlation', method="average", row_colors = row_colors, standard_scale=1, cmap='pink',figsize=(10,10))
sns_plot4.savefig("heatmap11-1.png")
plt.show()

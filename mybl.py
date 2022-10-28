import pandas as pd
import os

# Function of rbh in reeciprocal blast results

def my_blast(path,file_name):


    # Load reciprocal results

    os.chdir(os.path.expanduser(path))
    df_reci= pd.read_csv(file_name, sep=' ', header=None)
    #print(df_reci)


    # Add column names

    df_reci.columns=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'slen','qlen']
    df_reci.pident = df_reci.pident.astype(float)
    df_reci.bitscore = df_reci.bitscore.astype(float)
    #print(df_reci)


    # Add a col named pairs with format : qseqid sseqid

    df_reci["pairs"] = df_reci.apply(lambda row: ' '.join(sorted([row["qseqid"], row["sseqid"]])), axis = 1)
    #print(df_reci)


    # For every same element in pairs calculate the mean of...

    d1=df_reci.groupby('pairs', as_index=False)['bitscore'].mean()
    d2=df_reci.groupby('pairs', as_index=False)['pident'].mean()

    # Add those columns in the original df with desc or asc order

    metrics=pd.merge(d1, d2,   on="pairs")
    metrics['qseqid'] = metrics['pairs'].str.split(' ').str[1]
    metrics=metrics.sort_values(['qseqid','bitscore','pident'], ascending=[True,False,False])
    #print(metrics)


    # Keep the first of the results for every qseqid (best hit)
    metrics=metrics.drop_duplicates(subset=['qseqid'], keep='first')
    metrics = metrics.reset_index(drop=True)
    #print(metrics)


    # Now, we only keep these best results from the initial reciprocal df
    df_new = df_reci[df_reci.pairs.isin(metrics.pairs)]
    del df_new['pairs']
    #print(df_new)
    df_new.to_csv(file_name,header=None,index=False, sep=' ')

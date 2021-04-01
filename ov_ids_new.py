import pandas as pd
import re
import os

def ov_ids(path,org,all_ids,overl_file,clean_file):

    # Make a table with all existing ids from reference genome
    os.chdir(os.path.expanduser(path))
    counter=-1
    df_ids = pd.DataFrame(columns=['ID'])
    with open(all_ids) as f:
        for i in f:
            counter+=1
            name = re.findall('>+(\w+.\d+)',i)[0]
            chromosome = re.findall('%s:+(\w+):'% org,i)[0]
            start = re.findall('%s:+\w+:(\d+):' % org,i)[0]
            end = re.findall('%s:+\w+:\d+:(\d+):'% org,i)[0]
            df_ids.loc[counter, 'ID'] = name
            df_ids.loc[counter, 'Chrom'] = chromosome
            df_ids.loc[counter, 'Start'] = start
            df_ids.loc[counter, 'End'] = end


    df_ids.Chrom = df_ids.Chrom.astype(str)
    df_ids.Start = df_ids.Start.astype(int)
    df_ids.End = df_ids.End.astype(int)
    df_ids = df_ids.sort_values(by=['Chrom','Start'])
    df_ids = df_ids.reset_index(drop=True)

    # Find overlapped ids
    df_names= df_ids.copy()

    len_old = len(df_ids)
    len_new = 0
    length = [len_old, len_new]
    while True:
        print('round')
        for ind in range(0,len(df_ids)):
            if ind < len(df_ids):
                lis = df_ids.index[(df_ids.loc[ind,'Chrom'] == df_ids['Chrom']) & (df_ids.loc[ind, 'Start'] >= df_ids['Start']) & (df_ids['End'] >= df_ids.loc[ind, 'Start'])].tolist()
                if len(lis)>1:
                    leng=[]
                    for i in lis:
                        leng.append(df_ids.loc[i, 'End'] - df_ids.loc[i, 'Start'])
                    lis.remove(lis[leng.index(max(leng))])
                    if ind in lis:
                        df_ids=df_ids.drop(lis)
                        df_ids = df_ids.reset_index(drop=True)
                        ind-=1
                    else:
                        df_ids=df_ids.drop(lis)
                        df_ids = df_ids.reset_index(drop=True)
            else:
                break

        length[1] = len(df_ids)
        if length[0] == length[1]:
            break
        else:
            length =[length[1], length[0]]



    df_over=df_names[~df_names.ID.isin(df_ids.ID)]
    print(len(df_names))
    print(len(df_over))

    overl=df_over.ID.tolist()
    with open(overl_file, "w") as output:
            output.write(str(overl))

    df_ids.to_csv(clean_file, sep='\t')
    column_values = df_ids.ID.values.ravel()
    unique_values =  pd.unique(column_values)
    print(len(unique_values))

import os
import numpy as np
from collections import defaultdict
import re
import itertools
import pandas as pd
import ast

def order(path_rec,reciprocal, path_ref_clean, clean_ref, path_sub_clean, clean_sub, path_save_CGO, path_save_nCGO):
    os.chdir(os.path.expanduser(path_rec))



    # Make dcitionary with keys the human proteins and values the hits of macaca proteins.
    data = np.loadtxt(reciprocal, dtype=str)
    data2=np.sort(data[:,0:2],axis=1)
    #print(data2[0,:])
    Diction = defaultdict(list)
    for element in data2:
        Diction[element[1]].append(element[0])

    new_Dic = {a:list(set(b)) for a, b in Diction.items()}

    pr_hum=[]
    pr_mac=[]
    for key, value in new_Dic.items():
        pr_hum.append(key)
        pr_mac.append(value)
    print(len(pr_mac))
    print(pr_hum[0:2])
    print(pr_mac[0:2])


    os.chdir(os.path.expanduser(path_ref_clean))
    df_human = pd.read_csv(clean_ref, sep="\t")
    print(df_human)

    os.chdir(os.path.expanduser(path_sub_clean))
    df_macaca = pd.read_csv(clean_sub, sep="\t")


    # It's time to find which pairs are CGO and nCGO

    cCGO=-1
    cnCGO=-1
    df_CGO= pd.DataFrame(columns=['Prot_human','Prot_macaca'])
    df_nCGO= pd.DataFrame(columns=['Prot_human','Prot_macaca'])
    seen=[]
    counter=range(len(pr_mac))
    for i in counter:
        if i%1000 == 0:
            print('1000')
        for p in pr_mac[i]:
            pr_h=pr_hum[i]
            pr_m=p
            if [pr_h,pr_m] in seen:
                break
            seen.append([pr_h,pr_m])
            ind_h=df_human[df_human['ID']==pr_h].index.values[0]
            ind_m=df_macaca[df_macaca['ID']==pr_m].index.values[0]
            if ind_h==0 and ind_m==0:
                listA1 = [i for i, el in enumerate(pr_hum) if df_human.loc[ind_h+1,'ID'] in el]
                listB1 = [i for i, el in enumerate(pr_mac) if df_macaca.loc[ind_m+1,'ID'] in el]
                if any(x in listA1 for x in listB1):
                       cCGO+=1
                       df_CGO.loc[cCGO, 'Prot_human'] = pr_h
                       df_CGO.loc[cCGO, 'Prot_macaca'] = pr_m
                else:
                    cnCGO+=1
                    df_nCGO.loc[cnCGO, 'Prot_human'] = pr_h
                    df_nCGO.loc[cnCGO, 'Prot_macaca'] = pr_m
            elif ind_h==0 and ind_m!=0 and ind_m!=len(df_macaca)-1:
                listA1 = [i for i, el in enumerate(pr_hum) if df_human.loc[ind_h+1,'ID'] in el]
                listB_1 = [i for i, el in enumerate(pr_mac) if df_macaca.loc[ind_m-1,'ID'] in el]
                listB1 = [i for i, el in enumerate(pr_mac) if df_macaca.loc[ind_m+1,'ID'] in el]
                if (any(x in listA1 for x in listB_1) or any(x in listA1 for x in listB1)) :
                       cCGO+=1
                       df_CGO.loc[cCGO, 'Prot_human'] = pr_h
                       df_CGO.loc[cCGO, 'Prot_macaca'] = pr_m
                else:
                    cnCGO+=1
                    df_nCGO.loc[cnCGO, 'Prot_human'] = pr_h
                    df_nCGO.loc[cnCGO, 'Prot_macaca'] = pr_m
            elif ind_h==0 and ind_m==len(df_macaca)-1:
                listA1 = [i for i, el in enumerate(pr_hum) if df_human.loc[ind_h+1,'ID'] in el]
                listB_1 = [i for i, el in enumerate(pr_mac) if df_macaca.loc[ind_m-1,'ID'] in el]
                if any(x in listA1 for x in listB_1):
                    cCGO+=1
                    df_CGO.loc[cCGO, 'Prot_human'] = pr_h
                    df_CGO.loc[cCGO, 'Prot_macaca'] = pr_m
                else:
                    cnCGO+=1
                    df_nCGO.loc[cnCGO, 'Prot_human'] = pr_h
                    df_nCGO.loc[cnCGO, 'Prot_macaca'] = pr_m
            elif ind_h!=0 and ind_h!=len(df_human)-1 and ind_m==0:
                listA_1 =[i for i, el in enumerate(pr_hum) if df_human.loc[ind_h-1,'ID'] in el]
                listA1 = [i for i, el in enumerate(pr_hum) if df_human.loc[ind_h+1,'ID'] in el]
                listB1 = [i for i, el in enumerate(pr_mac) if df_macaca.loc[ind_m+1,'ID'] in el]
                if (any(x in listA1 for x in listB1) or any(x in listA_1 for x in listB1)):
                       cCGO+=1
                       df_CGO.loc[cCGO, 'Prot_human'] = pr_h
                       df_CGO.loc[cCGO, 'Prot_macaca'] = pr_m
                else:
                    cnCGO+=1
                    df_nCGO.loc[cnCGO, 'Prot_human'] = pr_h
                    df_nCGO.loc[cnCGO, 'Prot_macaca'] = pr_m
            elif ind_h!=0 and ind_h!=len(df_human)-1 and ind_m!=0 and ind_m!=len(df_macaca)-1:
                listA_1 =[i for i, el in enumerate(pr_hum) if df_human.loc[ind_h-1,'ID'] in el]
                listA1 = [i for i, el in enumerate(pr_hum) if df_human.loc[ind_h+1,'ID'] in el]
                listB_1 = [i for i, el in enumerate(pr_mac) if df_macaca.loc[ind_m-1,'ID'] in el]
                listB1 = [i for i, el in enumerate(pr_mac) if df_macaca.loc[ind_m+1,'ID'] in el]

                if (any(x in listA_1 for x in listB_1) or any(x in listA_1 for x in listB1) or any(x in listA1 for x in listB_1) or any(x in listA1 for x in listB1)):
                    cCGO+=1
                    df_CGO.loc[cCGO, 'Prot_human'] = pr_h
                    df_CGO.loc[cCGO, 'Prot_macaca'] = pr_m
                else:
                    cnCGO+=1
                    df_nCGO.loc[cnCGO, 'Prot_human'] = pr_h
                    df_nCGO.loc[cnCGO, 'Prot_macaca'] = pr_m
            elif ind_h!=0 and ind_h!=len(df_human)-1 and ind_m==len(df_macaca)-1:
                listA_1 =[i for i, el in enumerate(pr_hum) if df_human.loc[ind_h-1,'ID'] in el]
                listA1 = [i for i, el in enumerate(pr_hum) if df_human.loc[ind_h+1,'ID'] in el]
                listB_1 = [i for i, el in enumerate(pr_mac) if df_macaca.loc[ind_m-1,'ID'] in el]
                if (any(x in listA1 for x in listB_1) or any(x in listA_1 for x in listB_1)):
                       cCGO+=1
                       df_CGO.loc[cCGO, 'Prot_human'] = pr_h
                       df_CGO.loc[cCGO, 'Prot_macaca'] = pr_m
                else:
                    cnCGO+=1
                    df_nCGO.loc[cnCGO, 'Prot_human'] = pr_h
                    df_nCGO.loc[cnCGO, 'Prot_macaca'] = pr_m
            elif ind_h==len(df_human)-1 and ind_m==0:
                listA_1 = [i for i, el in enumerate(pr_hum) if df_human.loc[ind_h-1,'ID'] in el]
                listB1 = [i for i, el in enumerate(pr_mac) if df_macaca.loc[ind_m+1,'ID'] in el]
                if any(x in listA_1 for x in listB1) :
                       cCGO+=1
                       df_CGO.loc[cCGO, 'Prot_human'] = pr_h
                       df_CGO.loc[cCGO, 'Prot_macaca'] = pr_m
                else:
                    cnCGO+=1
                    df_nCGO.loc[cnCGO, 'Prot_human'] = pr_h
                    df_nCGO.loc[cnCGO, 'Prot_macaca'] = pr_m
            elif ind_h==len(df_human)-1 and ind_m==len(df_macaca)-1:
                listA_1 = [i for i, el in enumerate(pr_hum) if df_human.loc[ind_h-1,'ID'] in el]
                listB_1 = [i for i, el in enumerate(pr_mac) if df_macaca.loc[ind_m-1,'ID'] in el]
                if any(x in listA_1 for x in listB_1):
                    cCGO+=1
                    df_CGO.loc[cCGO, 'Prot_human'] = pr_h
                    df_CGO.loc[cCGO, 'Prot_macaca'] = pr_m
                else:
                    cnCGO+=1
                    df_nCGO.loc[cnCGO, 'Prot_human'] = pr_h
                    df_nCGO.loc[cnCGO, 'Prot_macaca'] = pr_m
            elif ind_h==len(df_human)-1 and ind_m!=0 and ind_m!=len(df_macaca)-1:
                listA_1 = [i for i, el in enumerate(pr_hum) if df_human.loc[ind_h-1,'ID'] in el]
                listB_1 = [i for i, el in enumerate(pr_mac) if df_macaca.loc[ind_m-1,'ID'] in el]
                listB1 = [i for i, el in enumerate(pr_mac) if df_macaca.loc[ind_m+1,'ID'] in el]
                if (any(x in listA_1 for x in listB_1) or any(x in listA_1 for x in listB1)) :
                    cCGO+=1
                    df_CGO.loc[cCGO, 'Prot_human'] = pr_h
                    df_CGO.loc[cCGO, 'Prot_macaca'] = pr_m
                else:
                    cnCGO+=1
                    df_nCGO.loc[cnCGO, 'Prot_human'] = pr_h
                    df_nCGO.loc[cnCGO, 'Prot_macaca'] = pr_m


    df_CGO.to_csv(path_save_CGO, header=None, index=None, sep=' ', mode='a')
    df_nCGO.to_csv(path_save_nCGO, header=None, index=None, sep=' ', mode='a')

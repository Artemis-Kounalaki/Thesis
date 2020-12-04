# Sequel of order 
up_down=1
listhum=[]
listmac=[]
a1='ENSP00000397118.2'
b1='ENSMMUP00000025509.4'
listhum.append(a1)
listmac.append(b1)
index_h=df_human[df_human['ID']==a1].index.values[0]
index_m=df_macaca[df_macaca['ID']==b1].index.values[0]
print(index_m,index_h)
while True:
    if index_h!=0 and index_h!=len(df_human)-1 and index_m!=0 and index_m!=len(df_macaca)-1:
        listA_1 =[i for i, el in enumerate(pr_hum) if df_human.loc[index_h-up_down,'ID'] in el]
        listA1 = [i for i, el in enumerate(pr_hum) if df_human.loc[index_h+up_down,'ID'] in el]
        listB_1 = [i for i, el in enumerate(pr_mac) if df_macaca.loc[index_m-up_down,'ID'] in el]
        listB1 = [i for i, el in enumerate(pr_mac) if df_macaca.loc[index_m+up_down,'ID'] in el]
        if (any(x in listA_1 for x in listB_1) or any(x in listA_1 for x in listB1) or any(x in listA1 for x in listB_1) or any(x in listA1 for x in listB1)):
            listhum.insert(len(listhum),df_human.loc[index_h+up_down,'ID'])
            listhum.insert(0,df_human.loc[index_h-up_down,'ID'])
            listmac.insert(len(listmac),df_macaca.loc[index_m+up_down,'ID'])
            listmac.insert(0,df_macaca.loc[index_m-up_down,'ID'])
            up_down+=1
        else:
            break
    elif index_h==0 and index_m==0:
        listA1 = [i for i, el in enumerate(pr_hum) if df_human.loc[index_h+up_down,'ID'] in el]
        listB1 = [i for i, el in enumerate(pr_mac) if df_macaca.loc[index_m+up_down,'ID'] in el]
        if any(x in listA1 for x in listB1):
            listhum.insert(len(listhum),df_human.loc[index_h+up_down,'ID'])
            listmac.insert(len(listmac),df_macaca.loc[index_m+up_down,'ID'])
            up_down+=1
        else:
            break
    elif index_h==0 and index_m!=0 and index_m!=len(df_macaca)-1:
        listA1 = [i for i, el in enumerate(pr_hum) if df_human.loc[index_h+up_down,'ID'] in el]
        listB_1 = [i for i, el in enumerate(pr_mac) if df_macaca.loc[index_m-up_down,'ID'] in el]
        listB1 = [i for i, el in enumerate(pr_mac) if df_macaca.loc[index_m+up_down,'ID'] in el]
        if (any(x in listA1 for x in listB_1) or any(x in listA1 for x in listB1)) :
            listhum.insert(len(listhum),df_human.loc[index_h+up_down,'ID'])
            listmac.insert(0,df_macaca.loc[index_m-up_down,'ID'])
            listmac.insert(len(listmac),df_macaca.loc[index_m+up_down,'ID'])
            up_down+=1
        else:
            break
    elif index_h==0 and index_m==len(df_macaca)-1:
        listA1 = [i for i, el in enumerate(pr_hum) if df_human.loc[index_h+up_down,'ID'] in el]
        listB_1 = [i for i, el in enumerate(pr_mac) if df_macaca.loc[index_m-up_down,'ID'] in el]
        if any(x in listA1 for x in listB_1):
            listhum.insert(len(listhum),df_human.loc[index_h+up_down,'ID'])
            listmac.insert(0,df_macaca.loc[index_m-up_down,'ID'])
            up_down+=1
        else:
            break
    elif index_h!=0 and index_h!=len(df_human)-1 and index_m==0:
        listA_1 =[i for i, el in enumerate(pr_hum) if df_human.loc[index_h-up_down,'ID'] in el]
        listA1 = [i for i, el in enumerate(pr_hum) if df_human.loc[index_h+up_down,'ID'] in el]
        listB1 = [i for i, el in enumerate(pr_mac) if df_macaca.loc[index_m+up_down,'ID'] in el]
        if (any(x in listA1 for x in listB1) or any(x in listA_1 for x in listB1)):
            listhum.insert(0,df_human.loc[index_h-up_down,'ID'])
            listhum.insert(len(listhum),df_human.loc[index_h+up_down,'ID'])
            listmac.insert(len(listmac),df_macaca.loc[index_m+up_down,'ID'])
            up_down+=1
        else:
            break
    elif index_h==len(df_human)-1 and index_m==0:
        listA_1 = [i for i, el in enumerate(pr_hum) if df_human.loc[index_h-up_down,'ID'] in el]
        listB1 = [i for i, el in enumerate(pr_mac) if df_macaca.loc[index_m+up_down,'ID'] in el]
        if any(x in listA_1 for x in listB1) :
            listhum.insert(0,df_human.loc[index_h-up_down,'ID'])
            listmac.insert(len(listmac),df_macaca.loc[index_m+up_down,'ID'])
            up_down+=1
        else:
            break
    elif index_h==len(df_human)-1 and index_m==len(df_macaca)-1:
        listA_1 = [i for i, el in enumerate(pr_hum) if df_human.loc[index_h-up_down,'ID'] in el]
        listB_1 = [i for i, el in enumerate(pr_mac) if df_macaca.loc[index_m-up_down,'ID'] in el]
        if any(x in listA_1 for x in listB_1):
            listhum.insert(0,df_human.loc[index_h-up_down,'ID'])
            listmac.insert(0,df_macaca.loc[index_m-up_down,'ID'])
            up_down+=1
        else:
            break

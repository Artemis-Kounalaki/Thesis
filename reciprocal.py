from collections import defaultdict
import numpy as np
import os
import itertools

# Function that finds from blast results the reciprocal results.

def reciprocal(path1, blast_results1,blast_results2, save1, save2) :

    #Read blast results - make numpy array sorted

    os.chdir(os.path.expanduser(path1))

    data1 = np.loadtxt(blast_results1, dtype=str)
    data2=np.sort(data1[:,0:2],axis=1)
    data11 =np.loadtxt(blast_results2, dtype=str)
    data22=np.sort(data11[:,0:2],axis=1)

    #Merge the two protein columns.

    c=np.core.defchararray.add(data2[:,0],data2[:,1])
    data1=np.column_stack((c,data1))

    c1=np.core.defchararray.add(data22[:,0],data22[:,1])
    data11=np.column_stack((c1,data11))

    #Make a dictionary with keys be the sorted merged transcript names.

    dict2=defaultdict(list)
    for i in data1:
        dict2[i[0]].append(list(i[1:]))

    Dict1=defaultdict(list)
    for i in data11:
        Dict1[i[0]].append(list(i[1:]))



    #Keep a dictionary with keys both appear
    #in dict2 and Dict1

    keep1={k:v for k,v in dict2.items() if k in Dict1}
    unkeep1=list(itertools.chain.from_iterable(keep1.values()))
    keep2={k:v for k,v in Dict1.items() if k in dict2}
    unkeep2=list(itertools.chain.from_iterable(keep2.values()))


    #Make a txt file with the transformed blast results.
    with open(save1, 'w') as out:
        for item in unkeep1:
            out.write(str(item)+'\n')
    out.close()
    with open(save2, 'w') as out:
        for item in unkeep2:
            out.write(str(item)+'\n')
    out.close()

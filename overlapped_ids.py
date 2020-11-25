import os
import re


def overlap(path,ids,org,blast_results,save_file):
    os.chdir(os.path.expanduser(path))
    myis=[]
    with open(ids) as f:
        for i in f:
            name = re.findall('>+(\w+.\d+)',i)[0]
            chromosome = re.findall('%s:+(\w+):'% org,i)[0]
            start = re.findall('%s:+\w+:(\d+):' % org,i)[0]
            end = re.findall('%s:+\w+:\d+:(\d+):'% org,i)[0]
            myis.append([name,chromosome ,start, end])


    id1=[]
    id2=[]
    with open(blast_results) as file:
        for line in file:
            id1.append(line.split('\t')[0])
            id2.append(line.split('\t')[1])
    pairs=list(zip(id1,id2))


    # in which list id is included

    overlapped_pr=[]
    for pair in pairs:
        a= [i for i, el in enumerate(myis) if pair[0] in el][0]
        b= [i for i, el in enumerate(myis) if pair[1] in el][0]

        ch1=myis[a][1]
        st1=int(myis[a][2])
        end1=int(myis[a][3])
        range1=range(st1,end1)

        ch2=myis[b][1]
        st2=int(myis[b][2])
        end2=int(myis[b][3])
        range2=range(st2,end2)


        if ch1==ch2 and ((st1>st2 and st1<end2) or (st2>st1 and st2<end1)):
            if (end1-st1)>(end2-st2):
                overlapped_pr.append(pair[0])
            elif (end2-st2)>(end1-st1):
                overlapped_pr.append(pair[1])
            #overlapped_pairs.append([pair[0],pair[1]])

    with open(save_file, "w") as output:
        output.write(str(overlapped_pr))

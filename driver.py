from alon_prune import alon_prune
import numpy as np
import pandas as pd
import json
import sys
from collections import Counter

def map_pol(num_pol):
    if num_pol==1:
        return "positive"
    elif num_pol==0:
        return "negative"
    else:
        return "unknown"

def label_results(alon_num):
    df = pd.read_csv(alon_num)
    df.columns = ["edge", "FM"]
    df.edge = df.edge.apply(eval)
    tf = pd.read_csv("progress.csv")
    tf["edge"] = tf[["source", "target"]].apply(tuple, axis=1)
    tf = tf.set_index("edge")
    df = df.set_index("edge")
    done = tf.join(df)
    done = done.sort_values(by="FM", ascending=False)
    done=done.reset_index()[["source","target","Polarity","FM"]]
    done.to_csv("labeled_results.csv")
    return 0


def parse_JSON(JSON):
    with open (JSON,"r") as f:
        data = json.load(f)

    titles = data[list(data.keys())[0]]['titles']
    source = data[list(data.keys())[0]]['interaction'][1]
    target = data[list(data.keys())[0]]['interaction'][0]
    polarity = data[list(data.keys())[0]]['interaction'][4]


    edge_mapping = dict()
    i = 1
    for title in titles:
        edge_mapping[i] = title
        i = i + 1

    s = pd.DataFrame(source, columns=["source"])
    t = pd.DataFrame(target, columns=["target"])
    p = pd.DataFrame(polarity, columns=["Polarity"])
    p["Polarity"]=p["Polarity"].apply(map_pol)

    edge_list = s.join(t)
    edge_list = edge_list.join(p)

    edge_list["source"] = edge_list["source"].map(edge_mapping)
    edge_list["target"] = edge_list["target"].map(edge_mapping)

    edge_list.to_csv("init_.csv",index=False)

def write_JSON(json_file_path, enrichment):

    df = pd.read_csv(r"progress.csv")
    f = open(json_file_path)
    filename= json_file_path.replace(".json","")

    ##load in old json file##
    data = json.load(f)
    ke=data.keys()

    mydict={}
    titles=data[list(ke)[0]]['titles']
    i=1
    reverse={}
    for title in titles:
        mydict[title]=i
        i=i+1

    edge_set_prime=list()
    conf=list()

    total_edges=0
    trimmed=0

    for x in range(0,len(data[list(ke)[0]]['interaction'][1])):
        total_edges=total_edges+1
        mt=(data[list(ke)[0]]['interaction'][0][x],data[list(ke)[0]]['interaction'][1][x])
        edge_set_prime.append(mt)
        if(enrichment=="strict"):
            conf.append(0)
        else:
            conf.append(-1)

    for row in df.iterrows():
    
        if (mydict[row[1][1]],mydict[row[1][0]]) in edge_set_prime:
            tupl=(mydict[row[1][1]],mydict[row[1][0]])
            ind=edge_set_prime.index(tupl)
            conf[ind]=1
            trimmed=trimmed+1

    data[list(ke)[0]]['interaction'][4]=conf

    ##write out file##
    fullpath=filename+"_"+"ENRICHED_"+".json"
    with open(fullpath, "w") as outfile:
        json.dump(data, outfile)
    outfile.close()
    print(fullpath)

    with open(filename+"_results.txt","w") as f:
        res="Out of "+str(total_edges)+" edges there were "+str(total_edges-trimmed)+" edges removed for a final edgecount of "+str(trimmed)
        f.write(res)
    f.close()


def main(argv):
    json_file_path= argv[1]
    parse_JSON(json_file_path)
    p = alon_prune("init_.csv")
    print("Parsed JSON")
    print("Begin Pruning")
    type_ = argv[3]
    init_motif= p.getMotifs()
    init_motif.to_csv("intitial_motif_comp.csv")

    if(argv[2]=="large"):
        p.large_net_prune(.10)
        print("Pruning Complete")
        ##enrich network##
        done=p.label_graph()
        ##label graph with results
        done=done[["Num_Alon"]]
        finished_file="Alon_number.csv"
        done.to_csv(finished_file)
        label_results(finished_file)
        write_JSON(json_file_path,type_)
        
    else:
        p.prune()
        print("Pruning Complete")
        ##enrich network##
        done=p.label_graph()
        ##label graph with results
        done=done[["Num_Alon"]]
        finished_file="Alon_number.csv"
        done.to_csv(finished_file)
        label_results(finished_file)
        write_JSON(json_file_path,type_)

main(sys.argv)
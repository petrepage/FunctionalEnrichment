import copy
from tkinter import E
import seaborn as sns
import numpy as np
import pandas as pd
import networkx as nx
from package.j_motif_finder import *



class alon_prune:

    def __init__(self,first_net):
        
        self.first_net=first_net
        self.G = nx.from_pandas_edgelist(pd.read_csv(first_net),"source","target",edge_attr=True,create_using=nx.DiGraph())
        self.pruned=""
            
            
    def safe_remove(self,graph,sink,source,edge_1,edge_2):
    
        graph.remove_edge(edge_1,edge_2)
        for node in graph.nodes:
            in_d=graph.in_degree(node)
            out_d=graph.out_degree(node)
            if(in_d==0 and out_d==0):
                return False
            if(in_d==0):
                if(node not in source):
                    return False
            if(out_d==0):
                if(node not in sink):
                    return False
        return True

    
    def prune(self):
        w = pd.DataFrame(columns=["iter","mean"])
        
        iter_=0

        sink=set()
        source=set()
        for node in self.G.nodes:
            in_d=self.G.in_degree(node)
            out_d=self.G.out_degree(node)
            if(in_d==0):
                source.add(node)
            if(out_d==0):
                sink.add(node)
            
        
        mmc=0

        
        ##make copy of network deep copy##
        temp = copy.deepcopy(self.G)
        count=0
        cont=True

        #set cont to true#
        x=j_motif_finder(self.first_net)
        final_d=x.get_motifs()

        ##create num alon, then sort##

        ##make a list of removed edges##
        removed_edges=[]

        ##testing##
        max_mean_alon=0
        me=0
        first_iter=final_d.Num_Alon.mean()
        w.loc[len(w.index)] = [0,first_iter] 
        while(me>=max_mean_alon):
            iter_=iter_+1
            
            if(iter_%10==0):
                nx.to_pandas_edgelist(temp).to_csv("progress.csv",index=False)
                w.loc[len(w.index)] = [iter_,me] 
            
            max_mean_alon=me
    
            kemp = copy.deepcopy(temp)
            Qemp = copy.deepcopy(temp)
    
            cont=True
            ##bridge detection##
    
            r_edge = final_d.tail(1).index.values[0]
        
            
            
            ##old stuff
    
            cont=self.safe_remove(kemp,sink,source,r_edge[0],r_edge[1])
            
    
            for oe in self.G.out_edges(r_edge[1]):
                if((oe[1] in sink) and self.G.in_degree(oe[1]) < 2):
                    cont=False
                    
            Qemp = Qemp.to_undirected()
            Qemp.remove_edge(r_edge[0],r_edge[1])
            l=nx.connected_components(Qemp)
            if(len(list(l)))>1:
                cont=False
            

            ##check if removal creates a disjoint##

            if(cont):
                temp.remove_edge(r_edge[0],r_edge[1])
                removed_edges.append(r_edge)
    
            ##drop from df 

            final_d=final_d.drop(index=r_edge)
            count=count+1
    
            ###recalc motif for new net##
    
            ##send to pandas edgelist
            q=nx.to_pandas_edgelist(temp)
            x.set_net(q)
            mm=x.get_motifs()
            me=mm.Num_Alon.mean()


        nx.to_pandas_edgelist(x.get_net()).to_csv("progress.csv",index=False)
        w.to_csv("alon_iters.csv",index=False)
        self.pruned=x.get_net()
        return x.get_net()

    def prune_curve(self):

        iter_ = 0

        sink = set()
        source = set()
        for node in self.G.nodes:
            in_d = self.G.in_degree(node)
            out_d = self.G.out_degree(node)
            if (in_d == 0):
                source.add(node)
            if (out_d == 0):
                sink.add(node)

        mmc = 0
        ged_df = pd.DataFrame(columns=["num_alon_mean"])

        ##make copy of network deep copy##
        temp = copy.deepcopy(self.G)
        count = 0
        cont = True

        # set cont to true#
        x = j_motif_finder(self.first_net)
        final_d = x.get_motifs()

        ##create num alon, then sort##

        ##make a list of removed edges##
        removed_edges = []

        ##testing##
        max_mean_alon = 0
        me = 0
        edges=len(self.G.edges)
        i=0
        while (i<edges):
            iter_ = iter_ + 1

            if (iter_ % 10 == 0):
                nx.to_pandas_edgelist(temp).to_csv("progress.csv", index=False)


            max_mean_alon = me

            kemp = copy.deepcopy(temp)
            Qemp = copy.deepcopy(temp)

            cont = True
            ##bridge detection##

            r_edge = final_d.tail(1).index.values[0]

            ##old stuff

            cont = self.safe_remove(kemp, sink, source, r_edge[0], r_edge[1])

            for oe in self.G.out_edges(r_edge[1]):
                if ((oe[1] in sink) and self.G.in_degree(oe[1]) < 2):
                    cont = False

            Qemp = Qemp.to_undirected()
            Qemp.remove_edge(r_edge[0], r_edge[1])
            l = nx.connected_components(Qemp)
            if (len(list(l))) > 1:
                cont = False

            ##check if removal creates a disjoint##

            if (cont):
                temp.remove_edge(r_edge[0], r_edge[1])
                removed_edges.append(r_edge)

            ##drop from df

            final_d = final_d.drop(index=r_edge)
            count = count + 1

            ###recalc motif for new net##

            ##send to pandas edgelist
            q = nx.to_pandas_edgelist(temp)
            x.set_net(q)
            mm = x.get_motifs()
            me = mm.Num_Alon.mean()
            ged_df.loc[len(ged_df.index)] = [me]
            i=i+1

        nx.to_pandas_edgelist(x.get_net()).to_csv("progress.csv", index=False)
        return ged_df


    def label_graph(self):
        test= j_motif_finder("progress.csv")
        df=test.get_motifs()
        return df

    def large_net_prune(self,percent):
        w = pd.DataFrame(columns=["iter","mean"])

        ##refactor prune to:
            ##remove edges based on a percentage of edges from the network
            ##can create orphans, sources, and sinks 
        iter_=0

        ##init sink and sources lists##
        sink=set()
        source=set()
        for node in self.G.nodes:
            in_d=self.G.in_degree(node)
            out_d=self.G.out_degree(node)
            if(in_d==0):
                source.add(node)
            if(out_d==0):
                sink.add(node)
            
        
        mmc=0

        
        ##make copy of network deep copy##
        temp = copy.deepcopy(self.G)
        count=0
        cont=True

        #set cont to true#
        x=j_motif_finder(self.first_net)
        final_d=x.get_motifs()
        first_iter=final_d.Num_Alon.mean()
        w.loc[len(w.index)] = [0,first_iter] 

        ##create num alon, then sort##

        ##make a list of removed edges##
        removed_edges=[]

        ##testing##
        max_mean_alon=0
        me=0

        ##get bottom 5% of motifs then remove safely##

        ##start the loop here based on some condition
        while(me>=max_mean_alon):
            print("batch:", iter_)
            last_net= copy.deepcopy(temp)
            ln=nx.to_pandas_edgelist(last_net)

            iter_=iter_+1
            
        ## write to csv for sue later ## nx.to_pandas_edgelist(temp).to_csv("progress.csv",index=False)
            
            max_mean_alon=me
    
            cont=True
            

        ##edge selected for removal

            num_rows = round(len(final_d)*(percent))

        ##list of edges for removal
            r_edge = list(final_d.tail(num_rows).index)
        
        ##check if removal creates an orphan

            for e in r_edge:
                kemp = copy.deepcopy(temp)
                Qemp = copy.deepcopy(temp)
                cont=self.safe_remove(kemp,sink,source,e[0],e[1])
    
                for oe in self.G.out_edges(e[1]):
                    if((oe[1] in sink) and self.G.in_degree(oe[1]) < 2):
                        cont=False
                    
        ##checks if removal causes a disconnect
                Qemp = Qemp.to_undirected()
                Qemp.remove_edge(e[0],e[1])
                l=nx.connected_components(Qemp)
                if(len(list(l)))>1:
                    cont=False
            

        ##check if removal creates a disjoint##


        ##removes edges if passes all other tests
                if(cont):
                    temp.remove_edge(e[0],e[1])
                    removed_edges.append(e)
    
            ##drop from df 

                final_d=final_d.drop(index=e)
                count=count+1
    
            ###recalc motif for new net##
    
            ##send to pandas edgelist
            q=nx.to_pandas_edgelist(temp)
            x.set_net(q)
            mm=x.get_motifs()
            me=mm.Num_Alon.mean()
            w.loc[len(w.index)] = [iter_,me] 
            w.to_csv("alon_iters.csv",index=False)

        x.set_net(ln)
        w.to_csv("alon_iters.csv",index=False)
        nx.to_pandas_edgelist(x.get_net()).to_csv("progress.csv",index=False)
        self.pruned=x.get_net()
        return x.get_net()
    

    def getMotifs(self):
        x=j_motif_finder(self.first_net)
        final_d=x.get_motifs()
        return final_d



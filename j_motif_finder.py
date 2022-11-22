import copy
import seaborn as sns
import numpy as np
import pandas as pd
import networkx as nx
from networkx import fast_gnp_random_graph


class j_motif_finder:

    def __init__(self, first_file):

        self.first_file = first_file
        self.first_net = nx.from_pandas_edgelist(pd.read_csv(self.first_file), "source", "target", edge_attr=True,
                                                 create_using=nx.DiGraph())

    def have_bidirectional_relationship(self, G, node1, node2):
        return G.has_edge(node1, node2) and G.has_edge(node2, node1)

    def set_net(self, first):
        self.first_net = nx.from_pandas_edgelist(first, "source", "target", edge_attr=True, create_using=nx.DiGraph())

    def get_net(self):
        return self.first_net

    def set_net_raw(self, graph):
        self.first_net = graph

    def motif_manage(self, alpha, beta, gamma, typ, ty, master_dict):

        ##alpha, beta, gamma are three nodes##

        ##edges to deal with are alpha beta, alpha gamma, and beta gamma##

        ab = (alpha, beta)
        ag = (alpha, gamma)
        bg = (beta, gamma)
        edgelist = [ab, ag, bg]
        for edge in edgelist:
            members = [alpha, beta, gamma, ty]
            if members not in master_dict[edge]:
                master_dict[edge].append(members)
                typ[edge] = typ[edge] + 1

    def motif_manage_two(self, alpha, beta, typ):
        ab = (alpha, beta)
        ba = (beta, alpha)
        typ[ab] = typ[ab] + 1
        typ[ba] = typ[ba] + 1

    def motif_manage_3(self, alpha, beta, gamma, typ, ty, master_dict):
        ##alpha, beta, gamma are three nodes##
        ##edges to deal with are alpha beta, alpha gamma, and beta gamma##
        ab = (alpha, beta)
        bc = (beta, gamma)
        ca = (gamma, alpha)
        edgelist = [ab, bc, ca]
        for edge in edgelist:
            members = {alpha, beta, gamma, ty}
            if members not in master_dict[edge]:
                master_dict[edge].append(members)
                typ[edge] = typ[edge] + 1

    def motif_check(self, graph, alpha, beta, gamma, md):
        t = md["inc"]
        ty = "none"
        if ((graph[alpha][beta]['Polarity'] == 'positive') and (graph[alpha][gamma]['Polarity'] == 'positive') and (
                graph[beta][gamma]['Polarity'] == 'positive')):
            t = md["C1"]
            ty = "C1"
        elif ((graph[alpha][beta]['Polarity'] == 'negative') and (graph[alpha][gamma]['Polarity'] == 'negative') and (
                graph[beta][gamma]['Polarity'] == 'positive')):
            t = md["C2"]
            ty = "C2"
        elif ((graph[alpha][beta]['Polarity'] == 'positive') and (graph[alpha][gamma]['Polarity'] == 'negative') and (
                graph[beta][gamma]['Polarity'] == 'negative')):
            t = md["C3"]
            ty = "C3"
        elif ((graph[alpha][beta]['Polarity'] == 'negative') and (graph[alpha][gamma]['Polarity'] == 'positive') and (
                graph[beta][gamma]['Polarity'] == 'negative')):
            t = md["C4"]
            ty = "C4"
        elif ((graph[alpha][beta]['Polarity'] == 'positive') and (graph[alpha][gamma]['Polarity'] == 'positive') and (
                graph[beta][gamma]['Polarity'] == 'negative')):
            t = md["I1"]
            ty = "I1"
        elif ((graph[alpha][beta]['Polarity'] == 'negative') and (graph[alpha][gamma]['Polarity'] == 'negative') and (
                graph[beta][gamma]['Polarity'] == 'negative')):
            t = md["I2"]
            ty = "I2"
        elif ((graph[alpha][beta]['Polarity'] == 'positive') and (graph[alpha][gamma]['Polarity'] == 'negative') and (
                graph[beta][gamma]['Polarity'] == 'positive')):
            t = md["I3"]
            ty = "I3"
        elif ((graph[alpha][beta]['Polarity'] == 'negative') and (graph[alpha][gamma]['Polarity'] == 'positive') and (
                graph[beta][gamma]['Polarity'] == 'positive')):
            t = md["I4"]
            ty = "I4"
        return t, ty

    @staticmethod
    def calc_ged(te, tn):
        dnodes = 0

        for node in te.nodes:
            if (node not in tn.nodes):
                tn.add_node(node)
                dnodes = dnodes + 1
        for node in tn.nodes:
            if (node not in te.nodes):
                te.add_node(node)
                dnodes = dnodes + 1

        ##calc ged##

        temp_a = nx.adjacency_matrix(te, nodelist=sorted(te.nodes()))
        temp_a1 = temp_a.todense()
        f_a = nx.adjacency_matrix(tn, nodelist=sorted(tn.nodes()))
        f_a1 = f_a.todense()
        res = np.subtract(f_a1, temp_a1)
        res = np.absolute(res)
        ged = np.sum(res) + dnodes

        return ged

    def get_motifs(self):

        master_motf_dict = {}
        master_3_motif_dict = {}

        G = self.first_net

        for edge in G.edges:
            master_motf_dict[edge] = []

        ##Dict for each motif_type##
        ##key is edge, value is membership count##
        C1 = {}

        ##Init Keys##
        C2 = {}
        C3 = {}
        C4 = {}
        I1 = {}
        I2 = {}
        I3 = {}
        I4 = {}
        inc = {}
        tdp = {}
        tdn = {}
        Fb = {}

        inc_t = {}
        Ref_Count = {}
        for edges in G.edges:
            C1[edges] = 0
            C2[edges] = 0
            C3[edges] = 0
            C4[edges] = 0
            I1[edges] = 0
            I2[edges] = 0
            I3[edges] = 0
            I4[edges] = 0
            inc[edges] = 0
            tdp[edges] = 0
            tdn[edges] = 0
            Fb[edges] = 0
            inc_t[edges] = 0
        md = {"C1": C1, "C2": C2, "C3": C3, "C4": C4, "I1": I1, "I2": I2, "I3": I3, "I4": I4, "inc": inc, "tdp": tdp,
              "tdn": tdn, "Fb": Fb, "inc_t": inc_t}

        biconnections = set()
        for u, v in G.edges():
            if u > v:  # Avoid duplicates, such as (1, 2) and (2, 1)
                v, u = u, v
            if self.have_bidirectional_relationship(G, u, v):
                biconnections.add((u, v))

        for u, v in biconnections:
            upol = G[u][v]['Polarity']
            vpol = G[v][u]['Polarity']
            if (upol == "negative"):
                if (vpol == "positive"):
                    self.motif_manage_two(u, v, tdn)
            if (upol == "positive"):
                if (vpol == "positive"):
                    self.motif_manage_two(u, v, tdp)
                if (vpol == "negative"):
                    self.motif_manage_two(u, v, tdn)
            else:
                self.motif_manage_two(u, v, inc_t)

        ##for each node in graph
        for node in G.nodes:
            first_node_reachable = []
            ##if there are at least two outgoing edges##
            ##for each directly reachable node, add to a list
            for edge in G.out_edges(node):
                first_node_reachable.insert(len(first_node_reachable), edge[1])
                for n in first_node_reachable:
                    for e in G.out_edges(n):
                        ##we have reached ffl##!
                        if (e[1] in first_node_reachable):
                            ##get the type of motif based on polarity##
                            typ, ty = self.motif_check(G, node, n, e[1], md)
                            ##insert the node
                            self.motif_manage(node, n, e[1], typ, ty, master_motf_dict)

        ## 3 node feedback loop Detection##
        for A_node in G.nodes:
            A_list = []
            for A_edge in G.out_edges(A_node):
                A_list.insert(len(A_list),A_edge[1])
                for B_node in A_list:
                    B_list=[]
                    for B_edge in G.out_edges(B_node):
                        B_list.insert(len(B_list),B_edge[1])
                        for C_node in B_list:
                            C_list=[]
                            for C_edge in G.out_edges(C_node):
                                C_list.insert(len(C_list),C_edge[1])
                                for X_node in C_list:
                                    if(X_node==A_node):
                                        typ= Fb
                                        ty="Fb"
                                        self.motif_manage_3(A_node,B_node,C_node,typ,ty,master_motf_dict)

        if(Ref_Count):
            for edge in G.edges:
                Ref_Count[edge] = G[edge[0]][edge[1]]["RefCount"]

        C1_df = pd.DataFrame.from_dict(C1, orient='index', columns=["C1"])
        C2_df = pd.DataFrame.from_dict(C2, orient='index', columns=["C2"])
        C3_df = pd.DataFrame.from_dict(C3, orient='index', columns=["C3"])
        C4_df = pd.DataFrame.from_dict(C4, orient='index', columns=["C4"])

        I1_df = pd.DataFrame.from_dict(I1, orient='index', columns=["I1"])
        I2_df = pd.DataFrame.from_dict(I2, orient='index', columns=["I2"])
        I3_df = pd.DataFrame.from_dict(I3, orient='index', columns=["I3"])
        I4_df = pd.DataFrame.from_dict(I4, orient='index', columns=["I4"])

        inc_df = pd.DataFrame.from_dict(inc, orient='index', columns=["inc"])

        tdn_df = pd.DataFrame.from_dict(tdn, orient='index', columns=["tdn"])
        tdp_df = pd.DataFrame.from_dict(tdp, orient='index', columns=["tdp"])

        fb_df = pd.DataFrame.from_dict(Fb, orient='index', columns=["fb"])
        inc_t_df = pd.DataFrame.from_dict(inc_t, orient='index', columns=["inc_t"])

        if(Ref_Count):
            Refcount_df = pd.DataFrame.from_dict(Ref_Count, orient='index', columns=["Ref_Count"])

        final_df = C1_df.join(C2_df).join(C3_df).join(C4_df).join(I1_df).join(I2_df).join(I3_df).join(I4_df).join(
            inc_df).join(tdn_df).join(tdp_df).join(fb_df).join(inc_t_df)

        non_mem = final_df.loc[final_df[final_df.columns].eq(0).all(1)].index.tolist()
        final_df["Num_Alon"] = final_df.sum(axis=1)

        final_df["Alon_Membership"] = final_df.index

        final_df["Alon_Membership"] = final_df["Alon_Membership"].apply(lambda x: 0 if x in non_mem else 1)

        if(Ref_Count):

            final_df = final_df.join(Refcount_df)

        final_df = final_df.sort_values(by=["Num_Alon"], ascending=False)

        return final_df
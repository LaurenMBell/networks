import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt
import time
import argparse
import pickle

#predictr/predictdir
#direction of change 

def vote(G, neighbor, target):
    #up = 1, down = -1
    n_dir = G.nodes[neighbor]["dir"]
    e_dir = G[neighbor][target]["dir"]
    v = n_dir * e_dir #(up/up or down/down = 1 | down/up or up/down = -1)
    return v


#G = [l0, l1, l2, l3...ln] = UNION OF ALL LAYERS

def define_next_layer(G, ln, visited):
    lm = set() #Ln+1

    for node in ln:
        neighbors = set(G.neighbors(node))
        for neighbor in neighbors:
            if neighbor not in visited: 
                #SHOULD ONLY BE NEIGHBORS IN THE NEXT LAYER, NOT VISITED OR CURRENT LAYER
                lm.add(neighbor)
                
    return lm

def reverse_puc(G, ln, visited, i, f=None, thresh=0.2, first=False):
    #function to take each node in a layer and find directionality for it
    to_remove_nodes = set() #set to collect nodes to remove instead of during iteration
    #to_remove_edges = set()

    lm = define_next_layer(G, ln, visited)

    for node in lm:
        if f: f.write(f"NODE: {node}\n")
        if f: f.write(f"LAYER: {i+1}\n")

        up = 0
        up_e = []
        down = 0
        down_e = []
        n_lm = []

        for neighbor in list(G.neighbors(node)):
            if neighbor in ln:
                n_lm.append(neighbor)
                #every neighbor in ln should vote 
                n_vote = vote(G, neighbor, node)
                if n_vote == 1:
                    up += 1
                    up_e.append(neighbor)
                if n_vote == -1:
                    down += 1
                    down_e.append(neighbor)

        score = up - down
        if f: f.write(f"SCORE: {score}\n")
        if f: f.write(f"UP: {up}\n")
        if f: f.write(f"DOWN: {down}\n")
        if score == 0: #tie or no neighbors
            frustration = None
            to_remove_nodes.add(node)
            if f: f.write("score == 0, node to be removed\n\n")

            if first:
                for n in n_lm:
                    G.remove_edge(node, n)
            continue
        elif score < 0:
            dir = -1
            frustration = up/(up+down)

            for neighbor in up_e:
                G.remove_edge(node, neighbor)
                n_lm.remove(neighbor)
                if f: f.write(f"removing edge {node} - {neighbor}\n")

        elif score > 0:
            dir = 1
            frustration = down/(up+down)

            for neighbor in down_e:
                G.remove_edge(node, neighbor)
                n_lm.remove(neighbor)
                if f: f.write(f"removing edge {node} - {neighbor}\n")

        if frustration is not None:
            if frustration > thresh: #if frustration is > 0.2, remove node
                to_remove_nodes.add(node)

                if first: 
                    for n in n_lm:
                        G.remove_edge(node, n)
                if f: f.write(f"{frustration} > {thresh}, node to be removed\n\n")
            else:
                G.nodes[node]['dir'] = dir
                if f: f.write("node updated!\n\n")

    if not first: G.remove_nodes_from(to_remove_nodes)

    return lm - to_remove_nodes

def pls_cpx_rpuc(f):
    f.write(time.strftime("CURRENT TIME: %Y-%m-%d %H:%M:%S\n\n"))
    f.write("Started PLS-CPX!\n")
    pls_cpx = pd.read_csv("PLS/PLS-CPX_edges.csv")
    pls_cpx = pls_cpx[pls_cpx["Pooled FDR"] <= 0.1]

    pls = pd.read_csv("PLS/PLS_edges.csv")
    pls = pls[(pls["FDR"] <= 0.05) & (pls[["VECPAC p-values", "DSS p-values", "LPS p-values"]].max(axis=1) <= 0.2)]
    cpx_node_dir = pd.read_csv("network_nodes.csv")
    
    l0 = set(pls_cpx['cpx_gene']) #going to be the cpx nodes in pls-cpx
    l1 =  set(pls_cpx['pls_metabolite']) #going to be pls nodes in pls-cpx

    #node direction hash table (gene:dir)
    dirs = dict(zip(cpx_node_dir["ID"], cpx_node_dir["Mean Log2 Fold Change Direction (DSS)"]))

    #init dir for all edges and cpx nodes
    G = nx.Graph()

    #init L0 and L1 from the PLS-CPX file
    for i, r in pls_cpx.iterrows():
        G.add_edge(r["cpx_gene"], r["pls_metabolite"], dir=r["edge_dir"]) 

    #add node DOC from netwrok_nodes
    for node in G.nodes: 
        if node in l0:
            G.nodes[node]["dir"] = dirs.get(node, None)

    #init edges for rest of PLS
    for i, r in pls.iterrows(): 
        G.add_edge(r["Metabolite 1"], r["Metabolite 2"], dir=r["edge_dir"])
    # ==================== PERFORMING INVERSE PUC =====================================
    # reverse_puc(G, ln, visited, f=None, thresh=0.2, first=False) -> lm - to_remove
    
    # GET L1 FROM L0
    visited = set(l0)

    l1 = reverse_puc(G, l0, visited, 0, f, first=True)

    visited |= l1 #visited = visited | l1

    layers = [set(l0), set(l1)] #list of sets of nodes

    ln = l1
    i = 1
    while True:
        lm = reverse_puc(G, ln, visited, i, f, first=False) #START WiTH L1 AS PLS 'L0'
        if not lm:
            #end of network
            break

        layers.append(lm)
        visited |= lm
        ln = lm #ln+1

        i+=1


    f.write("FINAL NETWORK NODES:\n")

    to_rmv = []
    for node in G.nodes:
        if "dir" not in G.nodes[node]:
            to_rmv.append(node)
        
    for n in to_rmv:
        G.remove_node(n)

    for node in G.nodes:
        f.write(f"{node}: {G.nodes[node]}\n")

    p = open("PLS/id_to_symbol_map.pickle", 'rb')
    name_dict = pickle.load(p)
    p.close()

    #save edge table 
    edges = []
    pls_e = []
    pls_cpx_e = []
    for u, v, d in G.edges(data=True):
        edges.append({"n1": u,"n2": v,"edge_dir": d.get("dir", None)})
        if u in l0 or v in l0:
            pls_cpx_e.append({"n1": u,"n2": v,"edge_dir": d.get("dir", None)})

        if u not in l0 and v not in l0:
            pls_e.append({"n1": u,"n2": v,"edge_dir": d.get("dir", None)})
    for u, v, d in G.edges(data=True):
        edges.append({"n1": v,"n2": u,"edge_dir": d.get("dir", None)})
        if u in l0 or v in l0:
            pls_cpx_e.append({"n1": v,"n2": u,"edge_dir": d.get("dir", None)})
        if u not in l0 and v not in l0:
            pls_e.append({"n1": v,"n2": u,"edge_dir": d.get("dir", None)})  
    
    pls_cpx_e = pd.DataFrame(pls_cpx_e)
    sel = pls_cpx_e["n1"].str.startswith("ENS")
    pls_cpx_e[sel].replace("-P", "")
    pls_cpx_e.loc[sel, "n1"] = pls_cpx_e.loc[sel, "n1"].str.replace("-P", "")
    sel2 = pls_cpx_e["n2"].str.startswith("ENS")
    pls_cpx_e[sel2].replace("-P", "")
    pls_cpx_e.loc[sel2, "n2"] = pls_cpx_e.loc[sel2, "n2"].str.replace("-P", "")

    pls_cpx_e = pls_cpx_e.replace(name_dict)
    pls_cpx_e.to_csv("PLS/pls_cpx_l0_rpuc_edges.csv", index=False)
    edges = pd.DataFrame(edges).replace("-P", "").replace(name_dict)
    edges.to_csv("PLS/pls_cpx_rpuc_edges.csv", index=False)
    pls_e = pd.DataFrame(pls_e).replace("-P", "").replace(name_dict)
    pls_e.to_csv("PLS/pls_rpuc_edges.csv", index=False)
    
    
    #save node table 
    nodes = []
    for n, d in G.nodes(data=True):
        nodes.append({"node": n,"node_dir": G.nodes[n]["dir"]})
       # if n not in l0:
            #try:
                #nodes.append({"node": n,"node_dir": G.nodes[n]["dir"]})
           # except:
                #print(f"{n} - no dir\n")
    pd.DataFrame(nodes).to_csv("PLS/pls_rpuc_nodes.csv", index=False) 

    count = 0
    for layer in layers:
        f.write(f"LAYER {count}: \n{layer}\n\n\n")
        count+=1

    G.clear()

def feci_cpx_rpuc(f):
    f.write(time.strftime("CURRENT TIME: %Y-%m-%d %H:%M:%S\n\n"))
    f.write("Started FECI-CPX!\n")
    feci_cpx = pd.read_csv("FECI/FECI-CPX_edges.csv")
    print("d-Erythrotetrofuranose" in feci_cpx["feci_metabolite"] ,"- d-Erythrotetrofuranose in feci-pls, expected FALSE")
    feci_cpx = feci_cpx[feci_cpx["Pooled FDR"] <= .1]

    feci = pd.read_csv("FECI/FECI_edges.csv")
    #print(("d-Erythrotetrofuranose" in feci["Metabolite 1"]) or ("d-Erythrotetrofuranose" in feci["Metabolite 2"]), "- expected TRUE")

    feci = feci[(feci["FDR"] <= 0.05) & (feci[["VECPAC p-values", "DSS p-values", "LPS p-values"]].max(axis=1) <= 0.2)]
    m = feci[
        (feci["Metabolite 1"] == "d-Erythrotetrofuranose") |
        (feci["Metabolite 2"] == "d-Erythrotetrofuranose")]
    print(m)
    
    print(f"{"d-Erythrotetrofuranose"} edges: {len(m)}")
        

    cpx_node_dir = pd.read_csv("network_nodes.csv")
    
    l0 = set(feci_cpx['cpx_gene']) #going to be the cpx nodes in feci-cpx
    l1 =  set(feci_cpx['feci_metabolite']) #going to be feci nodes in feci-cpx
    print("d-Erythrotetrofuranose" in l1, "- d-Erythrotetrofuranose in l1, expected FALSE")


    #node direction hash table (gene:dir)
    dirs = dict(zip(cpx_node_dir["ID"], cpx_node_dir["Mean Log2 Fold Change Direction (DSS)"]))

    #init dir for all edges and cpx nodes
    G = nx.Graph()

    #init L0 and L1 from the FECI-CPX file
    for i, r in feci_cpx.iterrows():
        G.add_edge(r["cpx_gene"], r["feci_metabolite"], dir=int(r["edge_dir"])) 

    #add node DOC from netwrok_nodes
    for node in G.nodes: 
        if node in l0:
            G.nodes[node]["dir"] = dirs.get(node, None)

    #init edges for rest of feci
    for i, r in feci.iterrows(): 
        G.add_edge(r["Metabolite 1"], r["Metabolite 2"], dir=int(r["edge_dir"]))

    print("d-Erythrotetrofuranose" in list(G.nodes), "- d-Erythrotetrofuranose in G.nodes, expected TRUE")


    # ==================== PERFORMING INVERSE PUC =====================================
    # reverse_puc(G, ln, visited, f=None, thresh=0.2, first=False) -> lm - to_remove
    
    # GET L1 FROM L0
    visited = set(l0)

    l1 = reverse_puc(G, l0, visited, 0, f, first=True)

    visited |= l1 #visited = visited | l1

    layers = [set(l0), set(l1)] #list of sets of nodes

    ln = l1
    i = 1
    while True:
        lm = reverse_puc(G, ln, visited, i, f, first=False) #START WiTH L1 AS PLS 'L0'
        if not lm:
            #end of network
            break

        layers.append(lm)
        visited |= lm
        ln = lm #ln+1

        i+=1

    f.write("FINAL NETWORK NODES:\n")

    to_rmv = []
    for node in G.nodes:
        if "dir" not in G.nodes[node]:
            to_rmv.append(node)
        
    for n in to_rmv:
        G.remove_node(n)

    for node in G.nodes:
        f.write(f"{node}: {G.nodes[node]}\n")

    try:
        for n in G.neighbors("d-Erythrotetrofuranose"):
            print(n)
    except:
        print("d-Erythrotetrofuranose not in network")

    

    #save edge table 
    edges = []
    feci_e = []
    for u, v, d in G.edges(data=True):
        edges.append({"n1": u,"n2": v,"edge_dir": d.get("dir", None)})
        if u not in l0 and v not in l0:
            feci_e.append({"n1": u,"n2": v,"edge_dir": d.get("dir", None)})
    pd.DataFrame(edges).to_csv("FECI/feci_cpx_rpuc_edges.csv", index=False)
    pd.DataFrame(feci_e).to_csv("FECI/feci_rpuc_edges.csv", index=False)
    
    
    #save node table 
    nodes = []
    for n, d in G.nodes(data=True):
        if n not in l0:
            try:
                nodes.append({"node": n,"node_dir": G.nodes[n]["dir"]})
            except:
                print(f"{n} - no dir\n")
    pd.DataFrame(nodes).to_csv("FECI/feci_rpuc_nodes.csv", index=False) 

    count = 0
    for layer in layers:
        f.write(f"LAYER {count}: \n{layer}\n\n\n")
        count+=1

    G.clear()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-pls", action="store_true")
    parser.add_argument("-feci", action="store_true")
    args = parser.parse_args()

    if args.pls:
        f = open("PLS/pls_report.txt", 'w')
        pls_cpx_rpuc(f)
    
    if args.feci:
        m = open("FECI/feci_report.txt", 'w')
        feci_cpx_rpuc(m)

if __name__=="__main__":
    main()

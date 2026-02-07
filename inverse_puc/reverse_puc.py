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
    n_dir = int(n_dir)
    e_dir = int(e_dir)
    v = n_dir * e_dir #(up/up or down/down = 1 | down/up or up/down = -1)
    return v


#G = [l0, l1, l2, l3...ln] = UNION OF ALL LAYERS

def build_layers(G, layer_0, max_depth=None):
    base_layer = []
    for node in layer_0:
        if G.has_node(node) and G.degree(node) > 0:
            base_layer.append(node)

    if not base_layer:
        return []

    layers = [set(base_layer)]
    visited = set(base_layer)
    current_layer = set(base_layer)
    depth = 0

    while current_layer:
        if max_depth is not None and depth >= max_depth:
            break

        next_layer = set()
        for node in current_layer:
            for neighbor in G.neighbors(node):
                if neighbor not in visited:
                    next_layer.add(neighbor)

        if not next_layer:
            break

        layers.append(next_layer)
        visited.update(next_layer)
        current_layer = next_layer
        depth += 1

    return layers

def reverse_puc(G, curr_layer, prev_layer, depth, f=None, thresh=0.2, first=False):
    for node in list(curr_layer):
        if not G.has_node(node):
            continue

        if f:
            f.write(f"\nNODE- {node},\nLAYER-{depth}\n")

        pos_votes = 0
        neg_votes = 0
        pos_edges = []
        neg_edges = []

        for neighbor in list(G.neighbors(node)):
            if neighbor not in prev_layer:
                continue
            if "dir" not in G.nodes[neighbor]:
                continue

            vote_val = vote(G, neighbor, node)
            if f:
                f.write(f"neighbor vote {neighbor}: {'up' if vote_val == 1 else 'down'}\n")

            if vote_val > 0:
                pos_votes += 1
                pos_edges.append((neighbor, node))
            else:
                neg_votes += 1
                neg_edges.append((neighbor, node))

        total_votes = pos_votes + neg_votes
        if f:
            f.write(f"positive votes: {pos_votes}\n") 
            f.write(f"negative votes: {neg_votes}\n")

        if total_votes == 0:
            print("NO. NO NO")
            if f:
                f.write("No votes, node unchanged\n")
            continue

        if pos_votes > neg_votes:
            G.nodes[node]["dir"] = 1
            minority_edges = neg_edges
        elif neg_votes > pos_votes:
            G.nodes[node]["dir"] = -1
            minority_edges = pos_edges
        else:
            minority_edges = pos_edges + neg_edges

        frustration = len(minority_edges) / total_votes

        if frustration > thresh:
            if first:
                if "dir" in G.nodes[node]:
                    del G.nodes[node]["dir"]
                for neighbor in list(G.neighbors(node)):
                    if neighbor in prev_layer and G.has_edge(node, neighbor):
                        G.remove_edge(node, neighbor)
                        if f:
                            f.write(f"{node}-{neighbor} edge is removed\n")
                if f:
                    f.write(f"{frustration} > {thresh}, edges cleared for {node}\n")
            else:
                if G.has_node(node):
                    G.remove_node(node)
                    if f:
                        f.write(f"{node} node is removed ({frustration} > {thresh})\n")
        else:
            for a, b in minority_edges:
                if G.has_edge(a, b):
                    G.remove_edge(a, b)
                    if f:
                        f.write(f"{a} - {b} edge is removed\n")

        if f:
            f.write(f"frustration={frustration}\n")

# deletes edges in same level
def same_level_edges(G, curr_layer, d, f):
    layer_nodes = set(n for n in curr_layer if G.has_node(n))

    for u in list(layer_nodes):
        for v in list(G.neighbors(u)):
            if v not in layer_nodes:
                continue
            if "dir" not in G.nodes[u] or "dir" not in G.nodes[v]:
                continue

            node_prod = float(G.nodes[u]["dir"]) * float(G.nodes[v]["dir"])
            edge_dir = float(G[u][v]["dir"])

            if node_prod != edge_dir and G.has_edge(u, v):
                if f:
                    f.write(f"SAME EDGE REMOVAL {v}-{u} edge is removed in layer {d}\n")
                G.remove_edge(u, v)

    return G




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
    layer_zero = [node for node in l0 if G.has_node(node) and G.degree(node) > 0]
    layers = build_layers(G, layer_zero)

    d = 1
    while d < len(layers):
        kept = 0
        for idx, layer in enumerate(layers):
            if idx > 0:
                kept += len(layer)

        prev_layer = set(layers[d - 1])
        curr_layer = set(layers[d])

        reverse_puc(G, curr_layer, prev_layer, d, f)
        same_level_edges(G, curr_layer, d, f)

        layers = build_layers(G, layer_zero)
        d += 1


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

    edges = pd.DataFrame(edges)
    sel = edges["n1"].str.startswith("ENS")
    edges[sel].replace("-P", "")
    edges.loc[sel, "n1"] = edges.loc[sel, "n1"].str.replace("-P", "")
    sel2 = edges["n2"].str.startswith("ENS")
    edges[sel2].replace("-P", "")
    edges.loc[sel2, "n2"] = edges.loc[sel2, "n2"].str.replace("-P", "")
    edges = edges.replace(name_dict)
    edges = pd.DataFrame(edges).replace("-P", "").replace(name_dict)
    edges.to_csv("PLS/pls_cpx_rpuc_edges.csv", index=False)

    pls_e = pd.DataFrame(pls_e)
    sel = pls_e["n1"].str.startswith("ENS")
    pls_e[sel].replace("-P", "")
    pls_e.loc[sel, "n1"] = pls_e.loc[sel, "n1"].str.replace("-P", "")
    sel2 = pls_e["n2"].str.startswith("ENS")
    pls_e[sel2].replace("-P", "")
    pls_e.loc[sel2, "n2"] = pls_e.loc[sel2, "n2"].str.replace("-P", "")
    pls_e = pls_e.replace(name_dict)
    pls_e = pd.DataFrame(pls_e).replace("-P", "").replace(name_dict)
    pls_e.to_csv("PLS/pls_rpuc_edges.csv", index=False)
    
    
    #save node table 
    nodes = []
    for n, d in G.nodes(data=True):
        #n = n.replace("-P", "").replace(name_dict)

        nodes.append({"node": n,"node_dir": G.nodes[n]["dir"]})
        if n not in l0:
            try:
                nodes.append({"node": n,"node_dir": G.nodes[n]["dir"]})
            except:
                print(f"{n} - no dir\n")
    pd.DataFrame(nodes).to_csv("PLS/pls_rpuc_nodes.csv", index=False) 

    count = 0
    for layer in layers:
        f.write(f"LAYER {count}: \n{layer}\n\n\n")
        count+=1

    G.clear()
    print("done with pls-cpx")

def fec_cpx_rpuc(f):
    f.write(time.strftime("CURRENT TIME: %Y-%m-%d %H:%M:%S\n\n"))
    f.write("Started FEC-CPX!\n")
    fec_cpx = pd.read_csv("FEC/FEC-CPX_edges_ananya.csv")
    fec_cpx = fec_cpx[fec_cpx["Pooled FDR"] <= 0.1]

    fec_pls = pd.read_csv("FEC/FEC-PLS_edges.csv")

    fec = pd.read_csv("FEC/FEC-FEC_edges_ananya.csv")
    fec = fec[(fec["FDR"] <= 0.05) & (fec[["VECPAC p-values", "DSS p-values", "LPS p-values"]].max(axis=1) <= 0.2)]
    cpx_node_dir = pd.read_csv("network_nodes.csv")
    
    l0 = set(fec_cpx['cpx_gene']) #going to be the cpx nodes in fec-cpx
    l1 =  set(fec_cpx['fec_metabolite']) #going to be fec nodes in fec-cpx

    #node direction hash table (gene:dir)
    dirs = dict(zip(cpx_node_dir["ID"], cpx_node_dir["Mean Log2 Fold Change Direction (DSS)"]))

    #init dir for all edges and cpx nodes
    G = nx.Graph()

    #init L0 and L1 from the FEC-CPX file
    for i, r in fec_cpx.iterrows():
        G.add_edge(r["cpx_gene"], r["fec_metabolite"], dir=r["edge_dir"]) 

    #add node DOC from netwrok_nodes
    for node in G.nodes: 
        if node in l0:
            G.nodes[node]["dir"] = dirs.get(node, None)

    #init edges for rest of FEC
    for i, r in fec.iterrows(): 
        G.add_edge(r["Metabolite 1"], r["Metabolite 2"], dir=r["edge_dir"])
    # ==================== PERFORMING INVERSE PUC =====================================
    layer_zero = [node for node in l0 if G.has_node(node) and G.degree(node) > 0]
    layers = build_layers(G, layer_zero)

    i = 1
    while i < len(layers):
        prev_layer = set(layers[i - 1])
        curr_layer = set(layers[i])

        reverse_puc(G, curr_layer, prev_layer, i, f)
        same_level_edges(G, curr_layer, i, f)

        layers = build_layers(G, layer_zero)
        i += 1

    f.write("FINAL NETWORK NODES:\n")

    to_rmv = []
    for node in G.nodes:
        if "dir" not in G.nodes[node]:
            to_rmv.append(node)
        
    for n in to_rmv:
        G.remove_node(n)

    for node in G.nodes:
        f.write(f"{node}: {G.nodes[node]}\n")

    p = open("FEC/id_to_symbol_map.pickle", 'rb')
    name_dict = pickle.load(p)
    p.close()

    #save edge table
    edges = []
    fec_e = []
    fec_cpx_e = []
    for u, v, d in G.edges(data=True):
        edges.append({"n1": u,"n2": v,"edge_dir": d.get("dir", None)})
        if u in l0 or v in l0:
            fec_cpx_e.append({"n1": u,"n2": v,"edge_dir": d.get("dir", None)})
        if u not in l0 and v not in l0:
            fec_e.append({"n1": u,"n2": v,"edge_dir": d.get("dir", None)})
    for u, v, d in G.edges(data=True):
        edges.append({"n1": v,"n2": u,"edge_dir": d.get("dir", None)})
        if u in l0 or v in l0:
            fec_cpx_e.append({"n1": v,"n2": u,"edge_dir": d.get("dir", None)})
        if u not in l0 and v not in l0:
            fec_e.append({"n1": v,"n2": u,"edge_dir": d.get("dir", None)})  
    
    fec_cpx_e = pd.DataFrame(fec_cpx_e)
    sel = fec_cpx_e["n1"].str.startswith("ENS")
    fec_cpx_e[sel].replace("-P", "")
    fec_cpx_e.loc[sel, "n1"] = fec_cpx_e.loc[sel, "n1"].str.replace("-P", "")
    sel2 = fec_cpx_e["n2"].str.startswith("ENS")
    fec_cpx_e[sel2].replace("-P", "")
    fec_cpx_e.loc[sel2, "n2"] = fec_cpx_e.loc[sel2, "n2"].str.replace("-P", "")
    fec_cpx_e = fec_cpx_e.replace(name_dict)
    fec_cpx_e.to_csv("FEC/fec_cpx_l0_rpuc_edges.csv", index=False)

    edges = pd.DataFrame(edges)
    sel = edges["n1"].str.startswith("ENS")
    edges[sel].replace("-P", "")
    edges.loc[sel, "n1"] = edges.loc[sel, "n1"].str.replace("-P", "")
    sel2 = edges["n2"].str.startswith("ENS")
    edges[sel2].replace("-P", "")
    edges.loc[sel2, "n2"] = edges.loc[sel2, "n2"].str.replace("-P", "")
    edges = edges.replace(name_dict)
    edges = pd.DataFrame(edges).replace("-P", "").replace(name_dict)
    edges.to_csv("FEC/fec_cpx_rpuc_edges.csv", index=False)

    fec_e = pd.DataFrame(fec_e)
    sel = fec_e["n1"].str.startswith("ENS")
    fec_e[sel].replace("-P", "")
    fec_e.loc[sel, "n1"] = fec_e.loc[sel, "n1"].str.replace("-P", "")
    sel2 = fec_e["n2"].str.startswith("ENS")
    fec_e[sel2].replace("-P", "")
    fec_e.loc[sel2, "n2"] = fec_e.loc[sel2, "n2"].str.replace("-P", "")
    fec_e = fec_e.replace(name_dict)
    fec_e = pd.DataFrame(fec_e).replace("-P", "").replace(name_dict)
    fec_e.to_csv("FEC/fec_rpuc_edges.csv", index=False)

    #save node table 
    nodes = []
    for n, d in G.nodes(data=True):
        nodes.append({"node": n,"node_dir": G.nodes[n]["dir"]})
        if n not in l0:
            try:
                nodes.append({"node": n,"node_dir": G.nodes[n]["dir"]})
            except:
                print(f"{n} - no dir\n")
    pd.DataFrame(nodes).to_csv("FEC/fec_rpuc_nodes.csv", index=False) 

    count = 0
    for layer in layers:
        f.write(f"LAYER {count}: \n{layer}\n\n\n")
        count+=1

    G.clear()
    print("done with fec-cpx")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-pls", action="store_true")
    parser.add_argument("-fec", action="store_true")
    args = parser.parse_args()

    if args.pls:
        f = open("PLS/pls_report.txt", 'w')
        pls_cpx_rpuc(f)
    
    if args.fec:
        m = open("FEC/fec_report.txt", 'w')
        fec_cpx_rpuc(m)

if __name__=="__main__":
    main()

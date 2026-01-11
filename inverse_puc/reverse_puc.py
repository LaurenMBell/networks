import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt
import igraph as ig
import time


#CAN THIS BE REWORKED WITH IGRAPH??

"""
1) define source of truth, 
    
    L0 = CPX NODES IN PLS-CPX NETWORK
    L1 = PLS NODES IN PLS-CPX NETWORK

2) for every node in Ln, every node that is directly connect to it that 
    isn't already in Ln-1 makes up the set Ln+1
3) Starting with L1, for each node X in set Ln, every node connected to it 
    in Ln-1. 
4) for each of those, they should vote on directionality with a sum mechanism
5) edges to nodes that disagree with the majority are removed 
5) frustration = proportion of disagreements, if that proportion is >0.2, 
    remove the node
"""

def vote(G, neighbor, target):
    #up = 1, down = -1
    n_dir = G.nodes[neighbor]["dir"]
    e_dir = G[neighbor][target]["dir"]
    v = n_dir * e_dir #(up/up or down/down = 1 | down/up or up/down = -1)
    return v


#G = [l0, l1, l2, l3...ln] = UNION OF ALL LAYERS

def define_layers(G, l0, l1):
    #YOU COULD DO THIS WITHOUT PASSING IN L1!!

    #function that returns a list of sets of nodes[l0,l1,...ln] 
    layers = [l0, l1]
    visited = l0 | l1

    while True:
        current = layers[-1]
        next_layer = set()

        for node in current:
            neighbors = set(G.neighbors(node))
            for neighbor in neighbors:
                if neighbor not in visited:
                    next_layer.add(neighbor)

        if not next_layer: #end at the end of the graph
            break

        layers.append(next_layer)
        visited |= next_layer # visited = visited | next_layer

    return layers


def reverse_puc(G, ln, lm, f, thresh=0.2):
    #function to take each node in a layer and find directionality for it
    to_remove_nodes = set() #set to collect nodes to remove instead of during iteration
    #to_remove_edges = set()

    for node in ln:
        f.write(f"NODE: {node}\n")

        up = 0
        up_e = []
        down = 0
        down_e = []

        for neighbor in list(G.neighbors(node)):
            if neighbor in lm:
                #every neighbor in lm should vote 
                n_vote = vote(G, neighbor, node)
                if n_vote == 1:
                    up += 1
                    up_e.append(neighbor)
                if n_vote == -1:
                    down += 1
                    down_e.append(neighbor)

        score = up - down
        f.write(f"SCORE: {score}\n")
        if score == 0: #tie or no neighbors
            to_remove_nodes.add(node)
            f.write("score == 0, node to be removed\n\n")
            continue
        elif score < 0:
            dir = -1
            frustration = up/(up+down)

            for neighbor in up_e:
                G.remove_edge(node, neighbor)
                f.write(f"removing edge {node} - {neighbor}\n")

        elif score > 0:
            dir = 1
            frustration = down/(up+down)

            for neighbor in down_e:
                G.remove_edge(node, neighbor)
                f.write(f"removing edge {node} - {neighbor}\n")

        

        if frustration > thresh: #if frustration is > 0.2, remove node
            to_remove_nodes.add(node)
            f.write("f > thresh, node to be removed\n\n")
        else:
            G.nodes[node]['dir'] = dir
            f.write("node updated!\n\n")

        
    
    G.remove_nodes_from(to_remove_nodes)

    return ln - to_remove_nodes

def pls_cpx_rpuc(f):
    f.write(f"CURRENT TIME: {time.localtime}\n\n")
    f.write("Started PLS-CPX!\n")
    pls_cpx = pd.read_csv("PLS-CPX_edges.csv")
    pls_cpx = pls_cpx[pls_cpx["Consistent"] == True]

    pls = pd.read_csv("PLS_edges.csv")
    pls = pls[pls["Consistent?"] == True]
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
    #TO DO NEED TO DOUBLE CHECK POOLED R CALCULATIONS FOR PLS EDGES !!! DO NOT FORGET TO DO THIS
    for i, r in pls.iterrows(): 
        G.add_edge(r["Metabolite 1"], r["Metabolite 2"], dir=r["edge_dir"])

    #define the rest of the graph
    layers = define_layers(G, l0, l1)
    
    #perform reverse_puc for all layers
    for i in range(1, len(layers)):
        layers[i] = reverse_puc(G, layers[i], layers[i-1], f)

    #to figure out later, it looks terrifying rn
    """
    #show the graph at the end
    nx.draw(G, with_labels=True, font_weight='bold')
    name = input("name the graph: ")
    plt.savefig(f"{name}.png")
    plt.show() 
    """ 

    f.write("FINAL NETWORK NODES:\n")
    for node in G.nodes:
        f.write(f"{node}: {G.nodes[node]}\n")

    #save edge table 
    edges = []
    for u, v, d in G.edges(data=True):
        edges.append({"n1": u,"n2": v,"edge_dir": d.get("dir", None)})
    pd.DataFrame(edges).to_csv("pls_cpx_rpuc_edges.csv", index=False)
    
    #save node table 
    nodes = []
    for n, d in G.nodes(data=True):
        if n not in l0:
            nodes.append({"node": n,"node_dir": G.nodes[n]["dir"]})
    pd.DataFrame(nodes).to_csv("pls_rpuc_nodes.csv", index=False)
    

def build_graph(edges):
    #function to make a graph from input edge table, and define L0 
    G = nx.Graph() #clearly just do this later 
    return G

def main():
    f = open("report.txt", 'w')
    pls_cpx_rpuc(f)

if __name__=="__main__":
    main()

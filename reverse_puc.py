import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt

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

        if not next_layer:
            break

        layers.append(next_layer)
        visited |= next_layer

    return layers


def reverse_puc(G, ln, lm):
    #function to take each node in a layer and find directionality for it
    to_remove = set() #list to collect nodes to remove instead of during iteration

    for node in ln:
        neighbors = list(G.neighbors(node))
        up = 0
        down = 0

        for neighbor in neighbors:
            if neighbor in lm:
                #every neighbor in lm should vote 
                n_vote = vote(G, neighbor, node)
                if n_vote == 1:
                    up += 1
                if n_vote == -1:
                    down += 1

        score = up - down
        if score < 0:
            G.nodes[node]['dir'] = -1
            f = up/(up+down)
        elif score > 0:
            G.nodes[node]['dir'] = 1
            f = down/(up+down)
        else: #score = 0, equal disagreement
            to_remove.add(node)
            continue

        if f >0.2: #if frustration is > 0.2, remove node
            to_remove.add(node)

    return ln - to_remove

def pls_cpx_rpuc():
    print("started pls-cpx")
    pls_cpx = pd.read_csv("PLS-CPX_edges.csv")
    pls = pd.read_csv("PLS_edges.csv")
    
    l0 = set(pls_cpx['cpx_gene']) #going to be the cpx nodes in pls-cpx
    l1 =  set(pls_cpx['pls_metabolite']) #going to be pls nodes in pls-cpx

    #you need to init 'dir' for all edges and nodes
    G = nx.Graph()

    for i, r in pls_cpx.iterrows():
        G.add_edge(r["cpx_gene"], r["pls_metabolite"], dir=r["edge_dir"])
    for i, r in pls.iterrows():
        G.add_edge(r["n1"], r["n2"], dir=r["edge_dir"])

    for node in G.nodes:
        G.nodes[node]["dir"] = G.nodes[node].get("dir", 1)

    layers = define_layers(G, l0, l1)
    
    for i in range(1, len(layers)):
        layers[i] = reverse_puc(G, layers[i], layers[i-1])

    #show the graph at the end
    nx.draw(G, with_labels=True, font_weight='bold')
    name = input("name the graph: ")
    plt.savefig(f"{name}.png")
    plt.show()

def main():
    pls_cpx_rpuc()

if __name__=="__main__":
    main()

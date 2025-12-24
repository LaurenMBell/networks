import networkx as nx
import pandas as pd

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
    #function that returns a list of subnetwork objects [l0,l1,...ln] 
    layers = [l0, l1]

    #something iterative needs to happen here
    #G = nx.union(G, current_layer) 

    for layer in layers:
        if layer == l0:
            continue
        else: #starting with l1
            next_layer = []

            for node in layer:
                neighbors = list(G.neighbors(node))
                for neighbor in neighbors:
                    if neighbor not in l0:
                        next_layer.append(neighbor)

            layers.append(next_layer)



    return layers


def reverse_puc(G, ln, lm):
    #function to take each node in a layer and find directionality for it

    # MAKE EACH LAYER A SET THEN FIND THE NEIGHBORS OF EACH NODE BELONING TO THAT SET
    for node in ln.nodes:
        neighbors = list(G.neighbors(node))
        score = 0

        #think of a better solution later 
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
                
                score += n_vote

        if score < 0:
            G.nodes[node]['dir'] = -1
            f = up/(up+down)
        elif score > 0:
            G.nodes[node]['dir'] = 1
            f = down/(up+down)
        else: #score = 0, equal disagreement
            G.remove_node(node)

        if f >0.2: #if frustration is > 0.2, remove node
            G.remove_node(node)

    return ln

"""
for every node in a layer:
    find every connecting node in the previous layer (all edges in ln-1)

    make each node vote and keep a running tally 

    if the tally is 0, throw it out

    if the tally is >0, the node is positive

    if the tally is <0, the node is negative 

after this works, find frustration
"""
        
def main():
    l0 = nx.Graph() #going to be the cpx nodes in pls-cpx
    l1 = nx.Graph() #going to be pls nodes in pls-cpx
    G = nx.union(l0, l1)

    layers = define_layers(G, l0, l1) #layers = list of subnetwork objects of G

    n = 0
    for layer in layers:
        if n == 0:
            continue
        else:
            layer = reverse_puc(G, layer, layers[n-1])
            n+=1
    

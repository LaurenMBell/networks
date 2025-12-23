import networkx as nx
import pandas as pd

"""
1) define source of truth, L0
2) for every node in Ln, every node that is directly connect to it that 
    isn't already in Ln-1 makes up the set Ln+1
3) Starting with L1, for each node X in set Ln, every node connected to it 
    in Ln-1. 
4) for each of those, they should vote on directionality with a sum mechanism
5) edges to nodes that disagree with the majority are removed 
5) frustration = proportion of disagreements, if that proportion is >0.2, 
    remove the node
"""


G = nx.Graph()

def reverse_puc(G, l, V):






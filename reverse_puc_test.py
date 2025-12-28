import pytest
import networkx as nx
import reverse_puc

#test all node/edge vote combinations 

def test_vote_up_up():
    G = nx.Graph()
    G.add_edge("A", "B", dir=1)
    G.nodes["A"]["dir"] = 1

    assert reverse_puc.vote(G, "A", "B") == 1

def test_vote_up_down():
    G = nx.Graph()
    G.add_edge("A", "B", dir=-1)
    G.nodes["A"]["dir"] = 1
    assert reverse_puc.vote(G, "A", "B") == -1


def test_vote_down_up():
    G = nx.Graph()
    G.add_edge("A", "B", dir=1)
    G.nodes["A"]["dir"] = -1
    assert reverse_puc.vote(G, "A", "B") == -1


def test_vote_down_down():
    G = nx.Graph()
   
    G.add_edge("A", "B", dir=-1)
    G.nodes["A"]["dir"] = -1
    assert reverse_puc.vote(G, "A", "B") == 1

def test_define_layers_3():
    G = nx.Graph()

    G.add_edge("A", "B", dir = 1)
    G.add_edge("B", "C", dir = 1)
    

    assert reverse_puc.define_layers()

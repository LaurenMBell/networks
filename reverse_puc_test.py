import pytest
import networkx as nx
import reverse_puc

#func for making a graph for testing
def make_graph(votes):
    #votes is a tuple
    G = nx.Graph()
    for n, t, n_dir, e_dir in votes:
        G.add_edge(n, t, dir=e_dir)
        G.nodes[n]["dir"] = n_dir
        G.nodes[t]["dir"] = 1 
    return G


# test all node combos
#++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
#+++++++++++++++++++++++++++++++++++++++++++++++++


def test_define_layers_straight():
    #for a straight graph, are all layers defined correctly?
    G = nx.Graph()
    G.add_edges_from([
        ("A", "B"),
        ("B", "C"),
        ("C", "D"),
    ])

    l0 = {"A"}
    l1 = {"B"}

    layers = reverse_puc.define_layers(G, l0, l1)
    assert layers == [{"A"}, {"B"}, {"C"}, {"D"}]

def test_define_layers_y_shape():
    #are layers defined correctly as a graph branches?
    G = nx.Graph()
    G.add_edges_from([
        ("A", "B"),
        ("B", "C"),
        ("B", "D"),
    ])

    l0 = {"A"}
    l1 = {"B"}

    layers = reverse_puc.define_layers(G, l0, l1)
    assert layers[2] == {"C", "D"}

def test_define_layers_stopping():
    #does the algorithm stop when it's supposed to?
    G = nx.Graph()
    G.add_edge("A", "B")
    G.add_edge("A", "C")

    l0 = {"A"}
    l1 = {"B", "C"}

    layers = reverse_puc.define_layers(G, l0, l1)
    assert len(layers) == 2

def test_define_layers_disconnected_graph():
    #does the algorithm skip any disconnencted components?
    G = nx.Graph()
    G.add_edges_from([
        ("A", "B"),
        ("C", "D"),
    ])

    l0 = {"A"}
    l1 = {"B"}

    layers = reverse_puc.define_layers(G, l0, l1)
    assert layers == [{"A"}, {"B"}]


def test_reverse_puc_up():
    #does reverse_puc vote right for a pos node?
    G = make_graph([
        ("A", "X",  1,  1), 
        ("B", "X",  1,  1), 
        ("C", "X", -1, -1),
        ("D", "X", 1, -1),
        ("E", "X", 1, 1),
        ("F", "X", -1, -1) 
    ])

    ln = {"X"}
    lm = {"A", "B", "C", "D", "E", "F"}

    remaining = reverse_puc.reverse_puc(G, ln, lm)

    assert "X" in remaining
    assert G.nodes["X"]["dir"] == 1

def test_reverse_puc_down():
    #does reverse_puc vote right for a neg node?
    G = make_graph([
        ("A", "X",  -1,  1), 
        ("B", "X",  1,  -1), 
        ("C", "X", -1, -1),
        ("D", "X", 1, -1),
        ("E", "X", -1, 1),
        ("F", "X", 1, -1) 
    ])

    ln = {"X"}
    lm = {"A", "B", "C", "D", "E", "F"}

    remaining = reverse_puc.reverse_puc(G, ln, lm)

    assert "X" in remaining
    assert G.nodes["X"]["dir"] == -1

def test_reverse_puc_frusteration_at_thresh():
    #shouldn't remove a node with 0.2 frustration
    G = make_graph([
        ("A", "X",  1, -1),
        ("B", "X", -1,  1),
        ("C", "X",  1,  1),
        ("D", "X",  -1,  1), 
        ("E", "X",  -1,  1)
    ])

    ln = {"X"}
    lm = {"A", "B", "C", "D", "E"}

    remaining = reverse_puc.reverse_puc(G, ln, lm)
    assert "X" in remaining


def test_reverse_puc_tie():
    #if score = 0, node should be removed
    G = make_graph([
        ("A", "X",  -1, -1),
        ("B", "X", -1,  1),
        ("C", "X",  1,  1),
        ("D", "X",  -1,  1)
    ])

    ln = {"X"}
    lm = {"A", "B", "C", "D"}

    remaining = reverse_puc.reverse_puc(G, ln, lm)
    assert remaining == set()

def test_reverse_puc_high_frustration():
    #frust at 0.3 should be removed
    G = make_graph([
        ("A", "X",  1,  1),
        ("B", "X",  1,  1),
        ("C", "X",  1, -1)
    ])
    #frustration 1/3 > 0.2

    ln = {"X"}
    lm = {"A", "B", "C"}

    remaining = reverse_puc.reverse_puc(G, ln, lm)

    assert "X" not in remaining

def test_reverse_puc_ignores_non_lm():
    #reverse_puc should only be pulling from lm, and NO WHERE ELSE
    G = make_graph([
        ("A", "X", 1,  1),
        ("B", "X", 1, -1),
        ("C", "X", 1, -1),
    ])

    ln = {"X"}
    lm = {"A"} #C AND B AREN"T IN PREV LAYER

    remaining = reverse_puc.reverse_puc(G, ln, lm)
    assert "X" in remaining
    assert G.nodes["X"]["dir"] == 1

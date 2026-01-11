import pytest
import networkx as nx
import reverse_puc

#func for making a graph for testing
def make_graph(votes):
    #votes is a tuple (n1, n2, n1_dir, n2_dir, e_dir)
    G = nx.Graph()
    for n1, n2, n1_dir, n2_dir, e_dir in votes:
        G.add_edge(n1, n2, dir=e_dir)
        G.nodes[n1]["dir"] = n1_dir
        G.nodes[n2]["dir"] = n2_dir #overwritten by reverse_puc 
        # ^^ SET TO NONE FOR YET TO BE PUCCED
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
    #for a straight graph, are all layers defined correctly
    G = nx.Graph()
    G.add_edges_from([("A", "B"), ("B", "C"), ("C", "D")]) 
    ##D - C - B - A

    l0 = {"A"}
    l1 = {"B"}

    layers = reverse_puc.define_layers(G, l0, l1)
    assert layers == [{"A"}, {"B"}, {"C"}, {"D"}]

def test_define_layers_y_shape():
    #are layers defined correctly as a graph branches
    # D\
    #   B - A
    # C/
    G = nx.Graph()
    G.add_edges_from([("A", "B"), ("B", "C"), ("B", "D")])

    l0 = {"A"}
    l1 = {"B"}

    layers = reverse_puc.define_layers(G, l0, l1)
    assert layers[2] == {"C", "D"}

def test_define_layers_2():
    # B\
    #   A
    # C/

    G = nx.Graph()
    G.add_edge("A", "B")
    G.add_edge("A", "C")

    l0 = {"A"}
    l1 = {"B", "C"}

    layers = reverse_puc.define_layers(G, l0, l1)
    assert len(layers) == 2 #only l0 and l1


def test_define_layers_disconnected_graph():
    #does the algorithm skip any disconnencted components
    G = nx.Graph()
    G.add_edges_from([("A", "B"),("C", "D")])

    l0 = {"A"}
    l1 = {"B"}

    layers = reverse_puc.define_layers(G, l0, l1)
    assert layers == [{"A"}, {"B"}]

def test_define_layers_complex_layers():
    #does the algorithm skip any connencted components
    #     E - D\
    #       \   B - A. -> and edge from A to E
    # G - F - C/  

    # layers should be L0: A, L1: B, E, L2: C, D, L3: F, L4, G
    G = nx.Graph()
    G.add_edges_from([("A", "B"), ("B", "C"), ("C", "F"), ("F", "G"), 
                      ("C", "E"), ("D", "E"), ("B", "D"), ("A", "E")])

    l0 = {"A"}
    l1 = {"B", "E"}

    layers = reverse_puc.define_layers(G, l0, l1)
    
    assert len(layers)==5


def test_reverse_puc_up():
    #does reverse_puc vote right for a pos node (small graph)
    G = make_graph([
        ("A", "X",  1,  None, 1), #p
        ("B", "X",  1,  None, 1), #p
        ("C", "X", -1, None, -1),#p
        ("D", "X", 1, None, -1),#n
        ("E", "X", 1, None, 1),#p
        ("F", "X", -1, None, -1)]) #p

    ln = {"X"}
    lm = {"A", "B", "C", "D", "E", "F"}

    remaining = reverse_puc.reverse_puc(G, ln, lm)

    assert "X" in remaining
    assert G.nodes["X"]["dir"] == 1

def test_reverse_puc_down():
    #does reverse_puc vote right for a neg node
    G = make_graph([
        ("A", "X",  -1,  None, 1), 
        ("B", "X",  1, None, -1), 
        ("C", "X", -1, None, -1),
        ("D", "X", 1, None, -1),
        ("E", "X", -1, None, 1),
        ("F", "X", 1, None, -1) ])

    ln = {"X"}
    lm = {"A", "B", "C", "D", "E", "F"}

    remaining = reverse_puc.reverse_puc(G, ln, lm)

    assert "X" in remaining
    assert G.nodes["X"]["dir"] == -1

def test_reverse_puc_frusteration_at_thresh():
    #shouldn't remove a node with 0.2 frustration
    G = make_graph([
        ("A", "X",  1, None, -1),
        ("B", "X", -1, None, 1),
        ("C", "X",  1, None, 1),
        ("D", "X",  -1, None, 1), 
        ("E", "X",  -1, None, 1)])

    ln = {"X"}
    lm = {"A", "B", "C", "D", "E"}

    remaining = reverse_puc.reverse_puc(G, ln, lm)
    assert "X" in remaining


def test_reverse_puc_tie():
    #if score = 0, node should be removed
    G = make_graph([
        ("A", "X",  -1, None, -1),
        ("B", "X", -1, None, 1),
        ("C", "X",  1, None, 1),
        ("D", "X",  -1, None, 1)])

    ln = {"X"}
    lm = {"A", "B", "C", "D"}

    remaining = reverse_puc.reverse_puc(G, ln, lm)
    assert remaining == set()

def test_reverse_puc_high_frustration():
    #frust at 0.3 should be removed
    G = make_graph([
        ("A", "X",  1, None, 1),
        ("B", "X",  1, None, 1),
        ("C", "X",  1, None, -1)])
    #frustration 1/3 > 0.2

    ln = {"X"}
    lm = {"A", "B", "C"}

    remaining = reverse_puc.reverse_puc(G, ln, lm)

    assert "X" not in remaining

def test_reverse_puc_ignores_non_lm():
    #reverse_puc should only be pulling from lm, and NO WHERE ELSE
    G = make_graph([
        ("A", "X", 1, None, 1),
        ("B", "X", 1, None, -1),
        ("C", "X", 1, None, -1)])

    ln = {"X"}
    lm = {"A"} #C AND B AREN"T IN PREV LAYER

    remaining = reverse_puc.reverse_puc(G, ln, lm)
    assert "X" in remaining
    assert G.nodes["X"]["dir"] == 1

#threshold as param (=0.2) X
#bigger graph (3+ layers and more nodes)
#network building ucntion (for bigger graph) 

def test_full_1():
    #big graph, 4 layers and 16 nodes, raised frustration threshold (for small graph)
    #votes is a tuple (n1, n2, n1_dir, n2_dir, e_dir)
    G = make_graph([
        #L3 to L2
        ("A", "E", None, None, 1), 
        ("B", "E", None, None, -1), 
        ("C", "E", None, None, -1), 
        ("B", "F", None, None, -1), 
        ("C", "F", None, None, -1), 
        ("D", "F", None, None, -1), #TO BE REMOVED, D TO BE REMOVED
        ("B", "G", None, None, -1), 
        ("C", "G", None, None, -1), #TO BE REMOVED
        ("D", "G", None, None, -1),  #TO BE REMOVED, D TO BE REMOVED

        #L2 to L1
        ("E", "H", None, None, 1), 
        ("F", "H", None, None, 1), 
        ("E", "I", None, None, -1), 
        ("F", "I", None, None, -1), 
        ("G", "I", None, None, 1), 
        ("G", "J", None, None, -1), #TO BE REMOVED
        ("I", "J", None, None, 1), #TO DO: need to account for within layer removal of disagreement 
        ("F", "K", None, None, -1), #TO BE REMOVED
        ("G", "K", None, None, -1), 

        #L1 to L0
        ("H", "L", None, 1, 1), 
        ("J", "L", None, 1, -1), 
        ("M", "H", None, 1, 1), 
        ("I", "M", None, 1, -1), 
        ("K", "M", None, 1, 1), 
        ("I", "N", None, -1, 1), 
        ("K", "N", None, -1, -1), 
        ("J", "O", None, 1, -1), 
        ("K", "O", None, 1, -1), 
        ("K", "P", None , -1, -1)
        ])

    ln = {"H", "I", "J", "K"}
    lm = {"L", "M", "N", "O", "P"}

    layers = reverse_puc.define_layers(G, lm, ln)

    assert len(layers) == 4

    for i in range(1, len(layers)):
        layers[i] = reverse_puc.reverse_puc(G, layers[i], layers[i-1], thresh = 0.5)

    #final predictions
    assert G.nodes["A"]["dir"] == 1 
    assert G.nodes["B"]["dir"] == -1
    assert G.nodes["C"]["dir"] == -1
    assert "D" not in G #D should be removed due to a tie in votes
    assert G.nodes["E"]["dir"] == 1 
    assert G.nodes["F"]["dir"] == 1 
    assert G.nodes["G"]["dir"] == -1
    assert G.nodes["H"]["dir"] == 1 
    assert G.nodes["I"]["dir"] == -1
    assert G.nodes["J"]["dir"] == -1
    assert G.nodes["K"]["dir"] == 1 


def test_full_2_high_f_thresh():
    #second big graph, 33 nodes and 5 layers, f thresh =0.5
    G = make_graph([ 
        #L0 to L1
        ("AA", "V", 1, None, 1),
        ("AA", "W", 1, None, -1),
        ("AA", "X", 1, None, 1),
        ("AA", "Z", 1, None, 1),
        ("AB", "V", 1, None, 1),
        ("AB", "W", 1, None, -1),
        ("AB", "X", 1, None, 1),
        ("AB", "Y", 1, None, -1),
        ("AC", "W", -1, None, 1),
        ("AC", "X", -1, None, -1),
        ("AC", "Y", -1, None, -1),
        ("AC", "Z", -1, None, -1),
        ("AD", "X", 1, None, -1),
        ("AD", "Y", 1, None, 1),
        ("AD", "Z", 1, None, 1),
        ("AE", "X", -1, None, 1),
        ("AE", "Z", -1, None, -1),
        ("AF", "Y", 1, None, 1),
        ("AF", "Z", 1, None, 1),
        ("AG", "Z", -1, None, 1),

        #L1 to L2
        ("V", "N", None, None, 1),
        ("V", "O", None, None, 1),
        ("V", "P", None, None, -1),
        ("V", "Q", None, None, -1),
        ("V", "R", None, None, 1),
        ("W", "N", None, None, -1),
        ("W", "O", None, None, 1),
        ("W", "P", None, None, 1),
        ("W", "Q", None, None, 1),
        ("X", "P", None, None, 1),
        ("X", "Q", None, None, -1),
        ("X", "R", None, None, 1),
        ("X", "S", None, None, 1),
        ("Y", "Q", None, None, 1),
        ("Y", "R", None, None, -1),
        ("Y", "S", None, None, -1),
        ("Y", "T", None, None, 1),
        ("Y", "U", None, None, -1),
        ("Z", "Q", None, None, 1),
        ("Z", "R", None, None, -1),
        ("Z", "S", None, None, -1),
        ("Z", "T", None, None, 1),
        ("Z", "U", None, None, -1),

        #L2 to L3
        ("N", "H", None, None, 1),
        ("N", "I", None, None, 1),
        ("N", "J", None, None, 1),
        ("O", "H", None, None, 1),
        ("O", "I", None, None, -1),
        ("O", "J", None, None, -1),
        ("O", "K", None, None, 1),
        ("O", "L", None, None, 1),
        ("P", "I", None, None, -1),
        ("P", "J", None, None, -1),
        ("P", "K", None, None, 1),
        ("P", "L", None, None, -1),
        ("P", "M", None, None, 1),
        ("Q", "J", None, None, 1),
        ("Q", "K", None, None, 1),
        ("Q", "L", None, None, -1),
        ("Q", "M", None, None, -1),
        ("R", "K", None, None, 1),
        ("R", "L", None, None, -1),
        ("R", "M", None, None, -1),
        ("R", "K", None, None, 1),
        ("S", "L", None, None, -1),
        ("S", "M", None, None, 1),
        ("T", "L", None, None, -1),
        ("T", "M", None, None, -1),
        ("U", "M", None, None, 1),

        #L3 to L4
        ("H", "A", None, None, 1),
        ("H", "B", None, None, 1),
        ("H", "C", None, None, 1),
        ("I", "A", None, None, -1),
        ("I", "B", None, None, 1),
        ("I", "C", None, None, -1),
        ("I", "D", None, None, -1),
        ("J", "B", None, None, -1),
        ("J", "C", None, None, -1),
        ("J", "D", None, None, -1),
        ("J", "B", None, None, -1),
        ("K", "C", None, None, 1),
        ("K", "D", None, None, -1),
        ("K", "E", None, None, 1),
        ("L", "D", None, None, -1),
        ("L", "E", None, None, -1),
        ("L", "F", None, None, 1),
        ("L", "G", None, None, -1),
        ("M", "E", None, None, 1),
        ("M", "F", None, None, 1),
        ("M", "G", None, None, 1)]) #holy network
    
    l0 = {"AA", "AB", "AC", "AD", "AE", "AF", "AG"}
    l1 = {"V", "W", "X", "Y", "Z"}

    layers = reverse_puc.define_layers(G, l0, l1)
    
    assert(len(layers) == 5)

    for i in range(1, len(layers)):
        print(f"\nLAYER: {i}")
        layers[i] = reverse_puc.reverse_puc(G, layers[i], layers[i-1], thresh = 0.5)
        print(f"\nFINAL NODES IN LAYER {i}: {layers[i]}")

    assert "A" not in G #tie
    assert G.nodes["B"]["dir"] == 1 
    assert G.nodes["C"]["dir"] == -1 
    assert G.nodes["D"]["dir"] == -1 
    assert G.nodes["E"]["dir"] == -1 
    assert "F" not in G #tie
    assert G.nodes["G"]["dir"] == -1 
    assert G.nodes["H"]["dir"] == 1 
    assert G.nodes["I"]["dir"] == 1 
    assert G.nodes["J"]["dir"] == 1 
    assert G.nodes["K"]["dir"] == -1 
    assert G.nodes["L"]["dir"] == 1 
    assert G.nodes["M"]["dir"] == -1 
    assert G.nodes["N"]["dir"] == 1 
    assert "O" not in G #tie 
    assert G.nodes["P"]["dir"] == -1 
    assert G.nodes["Q"]["dir"] == -1 
    assert "R" not in G #tie
    assert G.nodes["S"]["dir"] == -1 
    assert G.nodes["T"]["dir"] == 1 
    assert G.nodes["U"]["dir"] == -1 
    assert G.nodes["V"]["dir"] == 1 
    assert G.nodes["W"]["dir"] == -1 
    assert G.nodes["X"]["dir"] == 1 
    assert G.nodes["Y"]["dir"] == 1 
    assert G.nodes["Z"]["dir"] == 1 

    #assert "AG" not in G #should be removed even though its in L0
    
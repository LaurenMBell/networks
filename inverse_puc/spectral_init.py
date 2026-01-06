import networkx as nx
import reverse_puc

def make_graph(votes):
    #votes is a tuple (n1, n2, n1_dir, n2_dir, e_dir)
    G = nx.Graph()
    for n1, n2, n1_dir, n2_dir, e_dir in votes:
        G.add_edge(n1, n2, dir=e_dir)
        G.nodes[n1]["dir"] = n1_dir
        G.nodes[n2]["dir"] = n2_dir #overwritten by reverse_puc 
        # ^^ SET TO NONE FOR YET TO BE PUCCED
    return G

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
    
for i in range(1, len(layers)):
    print(f"\nLAYER: {i}")
    layers[i] = reverse_puc.reverse_puc(G, layers[i], layers[i-1], thresh = 0.5)
    print(f"\nFINAL NODES IN LAYER {i}: {layers[i]}")



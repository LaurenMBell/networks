#!/usr/bin/env python

import pandas as pd
import networkx as nx
import chime as c

# Returns a dictionary: node name -> BiBC
def bibc(G, nodes_0, nodes_1, normalized):
	bibcs = {n: 0.0 for n in G}
	for s in nodes_0:
		for t in nodes_1:
			# betweenness centrality does not count the endpoints (v not in s,t)
			paths_st = [x for x in nx.all_shortest_paths(G, s, t) if len(x) > 2]
			n_paths = len(paths_st)
			for path in paths_st:
				for n in path[1:-1]:  # Exclude endpoints
					bibcs[n] += 1 / n_paths

	if normalized:
		possible_paths_t0 = (len(nodes_0) - 1) * len(nodes_1)
		possible_paths_t1 = len(nodes_0) * (len(nodes_1) - 1)
		possible_paths_other = len(nodes_0) * len(nodes_1)
		for n in bibcs:
			if n in nodes_0:
				bibcs[n] /= possible_paths_t0
			elif n in nodes_1:
				bibcs[n] /= possible_paths_t1
			else:
				bibcs[n] /= possible_paths_other

	return bibcs


def make_edges():
	network_edges = pd.read_excel("gut_brain_network_2026_01_07.xlsx", sheet_name="Edges")

	fec_edges = pd.read_csv("FECI_RPUC_edges.csv")
	fec_edges["n1"] = fec_edges["n1"].str.strip()
	fec_edges["n2"] = fec_edges["n2"].str.strip()
	#feci_edges.loc[~feci_edges["n1"].str.startswith("ENSMUS"), "n1"] += '-F'
	#feci_edges.loc[~feci_edges["n2"].str.startswith("ENSMUS"), "n2"] += '-F'
	#feci_edges = feci_edges[(feci_edges["FDR"] <= 0.05) & (feci_edges["max_p"] <= 0.2)]

	fec_pls_edges = pd.read_csv("FECI-PLS_edges.csv")
	fec_pls_edges = fec_pls_edges[fec_pls_edges["FDR"] <= 0.1]
	fec_pls_edges["n1"] = fec_pls_edges["n1"].str.strip() + "-F"
	fec_pls_edges["n2"] = fec_pls_edges["n2"].str.strip()+ "-F"
	#feci_pls_edges["n1"] = feci_pls_edges["n1"].str.replace("-P", "")
	#feci_pls_edges["n2"] = feci_pls_edges["n1"].str.replace("-F", "")

	pls_edges = pd.read_csv("PLS_RPUC_edges.csv")
	pls_edges["n1"] = pls_edges["n1"].str.strip()
	pls_edges["n2"] = pls_edges["n2"].str.strip()
	pls_edges.loc[~pls_edges["n1"].str.startswith("ENSMUS"), "n1"] += '-P'
	pls_edges.loc[~pls_edges["n2"].str.startswith("ENSMUS"), "n2"] += '-P'
	#pls_edges = pls_edges[(pls_edges["FDR"] <= 0.05) & (pls_edges["max_p"] <= 0.2)]

	pls_cpx_edges = pd.read_csv("PLS-CPX_edges.csv")
	pls_cpx_edges = pls_cpx_edges[pls_cpx_edges["FDR"] <= 0.1]
	#pls_cpx_edges["n1"] = pls_cpx_edges["n1"].str.replace("-P", "")
	pls_cpx_edges["n1"] = pls_cpx_edges["n1"].str.strip()
	pls_cpx_edges["n2"] = pls_cpx_edges["n2"].str.strip() + "-P"

	cpx_edges = network_edges[(network_edges["Node 1 Tissue"] == "Choroid Plexus") |
		(network_edges["Node 2 Tissue"] == "Choroid Plexus")]
	
	cpx_edges = cpx_edges.rename(columns={"Node 1 Name": "n1", "Node 2 Name": "n2"})

	cpx_ctx_edges = network_edges[((network_edges["Node 1 Tissue"] == "Cortex") &
		 						(network_edges["Node 2 Tissue"] == "Choroid Plexus")) |
								((network_edges["Node 2 Tissue"] == "Cortex") &
		 						(network_edges["Node 1 Tissue"] == "Choroid Plexus"))]
	cpx_ctx_edges = cpx_ctx_edges.rename(columns={"Node 1 Name": "n1", "Node 2 Name": "n2"})

	ctx_edges = network_edges[(network_edges["Node 1 Tissue"] == "Cortex") |
		(network_edges["Node 2 Tissue"] == "Cortex")]
	ctx_edges = ctx_edges.rename(columns={"Node 1 Name": "n1", "Node 2 Name": "n2"})

	#edges is a list of tuples (edge table, sourcce tag)
	edges = [(fec_pls_edges, "FEC-PLS"), (pls_cpx_edges, "PLS-CPX"), (cpx_ctx_edges, "CPX-CTX"),
		(fec_edges, "FEC"), (pls_edges, "PLS"), (cpx_edges, "CPX"), (ctx_edges, "CTX")]
	
	return edges


def add_node_source(G, node, source):
	#add source tag 
	if not G.has_node(node):
		G.add_node(node, sources={source})
	else:
		G.nodes[node].setdefault("sources", set()).add(source)


def make_graph(G, edges):
    prenetwork = pd.DataFrame(columns=["n1", "n2", "source"])
    rows = [] 
	
    for table, source in edges:
        for i, r in table.iterrows():
            u, v = r["n1"], r["n2"]
            if not G.has_edge(v, u):
                G.add_edge(u, v, source=source)
            rows.append({"n1": u, "n2": v, "source": source}) #infiles manually changed to have n1, n2

			#add source for both nodes
            add_node_source(G, u, source)
            add_node_source(G, v, source)
    
    prenetwork = pd.DataFrame(rows)
    prenetwork.to_csv("preBiBC_edge_table.csv", index=False)
    return G


def nodes_in_source(G, source):
	#find all nodes with a given source tag (find tissue subnetworks)
	return {n for n, d in G.nodes(data=True)
		if source in d.get("sources", set())}

# FECI - PLS - CPX - CTX
edges = make_edges()
G = nx.Graph()
G = make_graph(G, edges)

print("PRE NODES: ", G.number_of_nodes())
print("PRE EDGES: ", G.number_of_edges())

giant_comp = max(nx.connected_components(G), key=len)
print("GIANT COMPONENT SIZE: ", len(giant_comp))

#N IS GIANT COMPONENT OF NETWORK
N = G.subgraph(giant_comp).copy()
#nx.write_edgelist(N, "giant_component_edge_table.csv")

print("POST NODES: ", N.number_of_nodes())
print("POST EDGES: ", N.number_of_edges())

print("CTX: ", len(nodes_in_source(N, "CTX")))
print("CPX: ", len(nodes_in_source(N, "CPX")))
print("PLS: ", len(nodes_in_source(N, "PLS")))
print("FEC: ", len(nodes_in_source(N, "FEC")))

fec_nodes = nodes_in_source(N, "FEC")
pls_nodes = nodes_in_source(N, "PLS")
ctx_nodes = nodes_in_source(N, "CTX")

bibc_fec_ctx = bibc(N, fec_nodes,  ctx_nodes, False)
bibc_pls_ctx = bibc(N, pls_nodes, ctx_nodes, False)

fec = (pd.DataFrame.from_dict(bibc_fec_ctx, orient="index", columns=["BiBC"]).reset_index(names="node"))
fec.to_csv("FEC-CTX_BiBC.csv", index=False)

pls = (pd.DataFrame.from_dict(bibc_pls_ctx, orient="index", columns=["BiBC"]).reset_index(names="node"))
pls.to_csv("PLS-CTX_BiBC.csv", index=False)
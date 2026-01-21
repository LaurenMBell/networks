import pandas as pd
import networkx as nx

# Returns a dictionary: node name -> BiBC
def bibc(G, nodes_0, nodes_1, normalized):
	bibcs = {n: 0.0 for n in G}
	for s in nodes_0:
		for t in nodes_1:
			# betweenness centrality does not count the endpoints (v not in s,t)
			paths_st = [x for x in list(nx.all_shortest_paths(G, s, t)) if len(x) > 2]
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

	feci_edges = pd.read_csv("FECI_edges.csv")
	feci_edges = feci_edges[(feci_edges["FDR"] <= 0.05) & (feci_edges["max_p"] <= 0.2)]

	feci_pls_edges = pd.read_csv("FECI-PLS_edges.csv")
	feci_pls_edges = feci_pls_edges[feci_pls_edges["FDR"] <= 0.1]

	pls_edges = pd.read_csv("PLS_edges.csv")
	pls_edges = pls_edges[(pls_edges["FDR"] <= 0.05) & (pls_edges["max_p"] <= 0.2)]

	pls_cpx_edges = pd.read_csv("PLS-CPX_edges.csv")
	pls_cpx_edges = pls_cpx_edges[pls_cpx_edges["FDR"] <= 0.1]

	#im assuming I can pull from the network table without any filtering 
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
	edges = [(feci_pls_edges, "FECI-PLS"), (pls_cpx_edges, "PLS-CPX"), (cpx_ctx_edges, "CPX-CTX"),
		(feci_edges, "FECI"), (pls_edges, "PLS"), (cpx_edges, "CPX"), (ctx_edges, "CTX")]
	return edges


def add_node_source(G, node, source):
	if not G.has_node(node):
		G.add_node(node, sources={source})
	else:
		G.nodes[node].setdefault("sources", set()).add(source)


def make_graph(G, edges):
	for table, source in edges:
		for i, r in table.iterrows():
			u, v = r["n1"], r["n2"] #each table manually changed to have n1, n2
			G.add_edge(u, v, source=source) #add edge to G

			#add source for both nodes
			add_node_source(G, u, source)
			add_node_source(G, v, source)
	return G


def nodes_in_source(G, source):
	#find all nodes with a given source tag (find subnetworks)
	return {n for n, d in G.nodes(data=True)
		if source in d.get("sources", set())}

# FECI - PLS - CPX - CTX
edges = make_edges()
G = nx.Graph()
G = make_graph(G, edges)

print("NODES: ", G.number_of_nodes())
print("EDGES: ", G.number_of_edges())

feci_nodes = nodes_in_source(G, "FECI")
pls_nodes = nodes_in_source(G, "PLS")
ctx_nodes = nodes_in_source(G, "CTX")

bibc_feci_ctx = bibc(G, feci_nodes, ctx_nodes, False)
bibc_pls_ctx = bibc(G, pls_nodes, ctx_nodes, False)

feci = (pd.DataFrame.from_dict(bibc_feci_ctx, orient="index", columns=["BiBC"]).reset_index(names="node"))
feci.to_csv("FECI-CTX_BiBC.csv", index=False)

pls = (pd.DataFrame.from_dict(bibc_pls_ctx, orient="index", columns=["BiBC"]).reset_index(names="node"))
pls.to_csv("PLS-CTX_BiBC.csv", index=False)

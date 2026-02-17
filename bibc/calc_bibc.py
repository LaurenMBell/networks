#!/usr/bin/env python
import argparse
from collections import defaultdict
import csv
from pathlib import Path
import networkx as nx
import numpy as np
import chime

def import_nw(nw_path):
	G = nx.Graph()
	with open(nw_path) as nw_file:
		csv_reader = csv.reader(nw_file)
		for row_index, row in enumerate(csv_reader):
			if row_index == 0:
				try:
					n1_index = int(row.index("node1"))
					n2_index = int(row.index("node2"))
				except ValueError:
					raise Exception("The input table must contain columns titled \"node1\" and \"node2\"")
			else:
				n1 = row[n1_index]
				n2 = row[n2_index]
				G.add_edge(n1, n2)
	return G

def read_node_type_map(type_map_path):
	type_map = defaultdict(list)
	with open(type_map_path) as type_map_file:
		csv_reader = csv.reader(type_map_file)
		for row in csv_reader:
			(node_name, node_type) = row
			type_map[node_type].append(node_name)
	return type_map

def subgraph_intersection(subgraph1, subgraph2):
	return np.intersect1d(subgraph1, subgraph2)

def get_type_subgraph(graph, type_map, type):
	if type not in type_map:
		raise Exception(f"Type {type} not in type map")
	return subgraph_intersection(type_map[type], graph)

def connected_component_subgraphs(network):
	return [network.subgraph(component) for component in nx.connected_components(network)]

# Returns a dictionary: node name -> BiBC
def bibc(G, nodes_0, nodes_1, normalized):
	bibcs = {n: 0.0 for n in G}
	for s in nodes_0:
		for t in nodes_1:
			# betweenness centrality does not count the endpoints (v not in s,t)
			paths_st = [x for x in list(nx.all_shortest_paths(G,s,t)) if len(x) > 2]
			n_paths = len(paths_st)
			for path in paths_st:
				for n in path[1:-1]: # Exclude endpoints
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

def write_bibc(bibc_map, file_path):
	with open(file_path, "w") as bibc_file:
		writer = csv.writer(bibc_file)
		writer.writerow(["node", "BiBC"])
		for n, bibc in bibc_map.items():
			writer.writerow([n, bibc])

parser = argparse.ArgumentParser(description="Computes the BiBC of every node in a network, relative to two types of nodes (i.e. subnetworks).")
parser.add_argument("--network", required=True, help="The path of a CSV file listing the edges in the network. It must include columns titled \"node1\" and \"node2\" that indicate the endpoints of each edge.")
parser.add_argument("--type_map", required=True, help="The path of a CSV file mapping nodes to types/classes. It must not have a header. The first column must contain node names and the second column must contain the type of the corresponding node.")
parser.add_argument("--type1", required=True, help="The node type to use for one end of the BiBC calculation.")
parser.add_argument("--type2", required=True, help="The node type to use for the other end of the BiBC calculation.")
parser.add_argument("--normalized", action="store_true", default=False, help="Normalize BiBC to be between 0 and 1. A normalized BiBC of 1 means the node appears on all shortest paths between the two node types; a value of 0 means it appears on none of them.")
parser.add_argument("--output", required=True, help="The path of a CSV file that will be produced containing the names of nodes (the first column) and their BiBCs (the second column).")
args = parser.parse_args()

network_path = Path(args.network)
if not network_path.exists:
	raise Exception(f"Specified network file {network_path} does not exist.")
	
output_path = Path(args.output)

G = import_nw(network_path)
gc = max(connected_component_subgraphs(G), key=len)

type_map = read_node_type_map(args.type_map)
type1_subgraph = get_type_subgraph(gc, type_map, args.type1)
type2_subgraph = get_type_subgraph(gc, type_map, args.type2)

print("loaded all data in")
bibc_map = bibc(gc, type1_subgraph, type2_subgraph, args.normalized)
write_bibc(bibc_map, output_path)
chime.success()
chime.success()
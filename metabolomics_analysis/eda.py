#Feb 2026, initial EDA of the networks 
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx

nodes = pd.read_excel("gut_brain_network_2026_02_17.xlsx", sheet_name="Nodes")
edges = pd.read_excel("gut_brain_network_2026_02_17.xlsx", sheet_name="Edges")

G = nx.from_pandas_edgelist(edges, source="Node 1 Name", target="Node 2 Name", create_using=nx.Graph())
nx.set_node_attributes(G, nodes.set_index("Node")["Tissue"].to_dict(), "tissue")

plt.figure(figsize=(8, 6))
for tissue, node_ids in nodes.groupby("Tissue")["Node"]:
    degrees = pd.Series(dict(G.degree(node_ids)))
    counts = degrees.value_counts().sort_index()
    plt.loglog(counts.index, counts.values, marker="o", linestyle="", label=tissue)

plt.xlabel("Node degree (log scale)")
plt.ylabel("Count (log scale)")
plt.title("Tissue-specific degree distributions")
plt.legend()
plt.tight_layout()
plt.show()
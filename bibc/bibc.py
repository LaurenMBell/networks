import pandas as pd
import networkx as nx
import chime as c

############################
# BiBC (unchanged logically)
############################

def bibc(G, nodes_0, nodes_1, normalized=False):
    bibcs = {n: 0.0 for n in G.nodes}

    for s in nodes_0:
        for t in nodes_1:
            if s == t:
                continue
            if not nx.has_path(G, s, t):
                continue

            paths = [p for p in nx.all_shortest_paths(G, s, t) if len(p) > 2]
            if not paths:
                continue

            w = 1 / len(paths)
            for p in paths:
                for n in p[1:-1]:
                    bibcs[n] += w

    return bibcs


########################################
# Utility: allowed nodes from base layer
########################################

def nodes_from_table(df):
    return set(df["n1"]) | set(df["n2"])


########################################
# Load + VALIDATE edge tables
########################################

def make_edges():

    # --- Base intranetwork layers (authoritative node sources)
    feci_edges = pd.read_csv("FECI_RPUC_edges.csv")
    pls_edges  = pd.read_csv("PLS_RPUC_edges.csv")

    # --- Gut–brain atlas (CPX / CTX)
    network_edges = pd.read_excel(
        "gut_brain_network_2026_01_07.xlsx",
        sheet_name="Edges"
    )

    cpx_edges = network_edges[
        (network_edges["Node 1 Tissue"] == "Choroid Plexus") |
        (network_edges["Node 2 Tissue"] == "Choroid Plexus")
    ].rename(columns={"Node 1 Name": "n1", "Node 2 Name": "n2"})

    ctx_edges = network_edges[
        (network_edges["Node 1 Tissue"] == "Cortex") |
        (network_edges["Node 2 Tissue"] == "Cortex")
    ].rename(columns={"Node 1 Name": "n1", "Node 2 Name": "n2"})

    # --- Internetwork edges (potentially dangerous)
    feci_pls = pd.read_csv("FECI-PLS_edges.csv")
    pls_cpx  = pd.read_csv("PLS-CPX_edges.csv")

    # ============================
    # DEFINE VALID NODE UNIVERSES
    # ============================

    feci_nodes = nodes_from_table(feci_edges)
    pls_nodes  = nodes_from_table(pls_edges)
    cpx_nodes  = nodes_from_table(cpx_edges)
    ctx_nodes  = nodes_from_table(ctx_edges)

    # ============================
    # HARD FILTER INTERNETWORKS
    # ============================

    feci_pls = feci_pls[
        feci_pls["n1"].isin(feci_nodes | pls_nodes) &
        feci_pls["n2"].isin(feci_nodes | pls_nodes)
    ]

    pls_cpx = pls_cpx[
        pls_cpx["n1"].isin(pls_nodes | cpx_nodes) &
        pls_cpx["n2"].isin(pls_nodes | cpx_nodes)
    ]

    # ============================
    # RETURN CLEAN EDGE SET
    # ============================

    edges = [
        (feci_edges, "FECI"),
        (pls_edges, "PLS"),
        (cpx_edges, "CPX"),
        (ctx_edges, "CTX"),
        (feci_pls, "FECI-PLS"),
        (pls_cpx, "PLS-CPX"),
    ]

    return edges


########################################
# Graph construction (NO node creation)
########################################

def make_graph(edges):
    G = nx.Graph()

    # add base layers first
    for table, source in edges:
        for _, r in table.iterrows():
            u, v = r["n1"], r["n2"]
            G.add_edge(u, v, source=source)

    return G


########################################
# Connectivity enforcement
########################################

def assert_connected(G):
    comps = list(nx.connected_components(G))
    components = sorted(nx.connected_components(G), key=len)
    for comp in components[:-1]:
        print(len(comp), comp)
        
    if len(comps) > 1:
        sizes = sorted(len(c) for c in comps)
        raise RuntimeError(
            f"Graph disconnected: {len(comps)} components "
            f"(largest={max(sizes)}, smallest={min(sizes)})"
        )


########################################
# Node selection
########################################

def nodes_in_layer(G, layer):
    return {
        n for n, d in G.nodes(data=True)
        if any(layer in e["source"] for _, _, e in G.edges(n, data=True))
    }


########################################
# MAIN
########################################

edges = make_edges()
G = make_graph(edges)

print("NODES:", G.number_of_nodes())
print("EDGES:", G.number_of_edges())

assert_connected(G)

feci_nodes = nodes_in_layer(G, "FECI")
pls_nodes  = nodes_in_layer(G, "PLS")
ctx_nodes  = nodes_in_layer(G, "CTX")

bibc_feci_ctx = bibc(G, feci_nodes, ctx_nodes)
bibc_pls_ctx  = bibc(G, pls_nodes, ctx_nodes)

pd.DataFrame.from_dict(
    bibc_feci_ctx, orient="index", columns=["BiBC"]
).reset_index(names="node").to_csv(
    "FECI-CTX_BiBC.csv", index=False
)

pd.DataFrame.from_dict(
    bibc_pls_ctx, orient="index", columns=["BiBC"]
).reset_index(names="node").to_csv(
    "PLS-CTX_BiBC.csv", index=False
)

c.success()

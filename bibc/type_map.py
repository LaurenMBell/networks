import pandas as pd

"""
1) read in CTX-CPX, CPX-CPX
2) add CPX-PLS and match genes to CPX and mets to PLS
3) add PLS-PLS and ret pls-ctx edge table and type map 
4) add PLS-FEC and do puc check for spurious edges
5) add FEC-FEC and ret fec-ctx edge table/typemap

- how to keep track of nodes throughout: dict with node info? df of node + tissues
- HOW TO DO FEC-PLS CHECK
"""
import pandas as pd
import numpy as np
import chime

# {node: {"tissues": set(), "ids": set()}}
nodes = {} 
tissue_flag = {"Plasma": "PLS","Feces": "FEC",
    "Choroid Plexus": "CPX", "Cortex": "CTX",
    "Striatum": "STM"}
flags = {v: k for k, v in tissue_flag.items()} #reverse tissue_flag

def load_node(nodes, name, tissue, node_id=None):
    #puts nodes into global nodes dictionary w/ tissues and ID if in excel (only for ctx)
    entry = nodes.setdefault(name, {"tissues": set(), "ids": set()})
    entry["tissues"].add(tissue)
    if node_id is not None: #PLS AND FEC HAVE NO ID
        entry["ids"].add(str(node_id))

def unknowns(name, tissue, node_id=None):
    name = str(name).strip()
    if name == "UNKNOWN" and tissue == "Cortex" and node_id is not None:
        return f"UNKNOWN_{node_id}"
    return name


def label_node(name, tissue):
    tag = "__" + tissue_flag[tissue]
    return name if name.endswith(tag) else name + tag


def dedup(df):
    #admittedly I got help from chatgpt for this one
    if df.empty:
        return df
    df = df.copy()
    df["min"] = df[["node1","node2"]].min(axis=1)
    df["max"] = df[["node1","node2"]].max(axis=1)
    df = df.drop_duplicates(subset=["min","max"])
    return df.drop(columns=["min","max"])

def load_edges(df, nodes, strip_tag=False):
    edges = []
    for i, row in df.iterrows():
        t1 = row["Node 1 Tissue"]
        t2 = row["Node 2 Tissue"]

        n1_pre = row["Node 1 Name"]
        n2_pre = row["Node 2 Name"]

        if strip_tag:
            n1_pre = strip_tissue_tag(n1_pre)
            n2_pre = strip_tissue_tag(n2_pre)

        n1 = unknowns(n1_pre, t1)
        n2 = unknowns(n2_pre, t2)

        load_node(nodes, n1, t1)
        load_node(nodes, n2, t2)

        edges.append({
            "node1": label_node(n1, t1),
            "node2": label_node(n2, t2)
        })

    return pd.DataFrame(edges)



def load_edges_excel(df, tissue1, tissue2, nodes, use_ids=False):
    df = df.rename(columns={"Node 1 Name": "node1", "Node 2 Name": "node2"}).copy()

    mask = ((df["Node 1 Tissue"] == tissue1) & (df["Node 2 Tissue"] == tissue2)) | (
        (df["Node 1 Tissue"] == tissue2) & (df["Node 2 Tissue"] == tissue1))

    df = df.loc[mask, ["node1","node2","Node 1 ID","Node 2 ID","Node 1 Tissue","Node 2 Tissue"]]

    edges = []
    for i, row in df.iterrows():
        n1 = unknowns(row["node1"], row["Node 1 Tissue"],
                        row["Node 1 ID"] if use_ids else None)
        n2 = unknowns(row["node2"], row["Node 2 Tissue"],
                    row["Node 2 ID"] if use_ids else None)

        load_node(nodes, n1, row["Node 1 Tissue"], row["Node 1 ID"])
        load_node(nodes, n2, row["Node 2 Tissue"], row["Node 2 ID"])

        edges.append({
            "node1": label_node(n1, row["Node 1 Tissue"]),
            "node2": label_node(n2, row["Node 2 Tissue"])})

    return dedup(pd.DataFrame(edges))

############### PUC FOR FEC-PLS ###################################
def strip_tissue_tag(name):
    #fec-pls still has -P and -F
    name = str(name).strip()
    if name.endswith("-F") or name.endswith("-P"):
        return name[:-2]
    return name


def build_dir_map(df):
    return {strip_tissue_tag(r["node"]): r["node_dir"] for _, r in df.iterrows()}


def puc_fec_pls(df):
    fec_df = pd.read_csv("edges/FEC_nodes.csv")
    pls_df = pd.read_csv("edges/PLS_nodes.csv")

    fec_dir = build_dir_map(fec_df)
    pls_dir = build_dir_map(pls_df)

    keep = []
    for i, row in df.iterrows():
        pls = strip_tissue_tag(row["Node 1 Name"])
        fec = strip_tissue_tag(row["Node 2 Name"])
        edge_dir = row["edge_dir"]

        if edge_dir is None:
            continue
        if pls not in pls_dir or fec not in fec_dir:
            continue
        if edge_dir * fec_dir[fec] == pls_dir[pls]:
            keep.append(row)

    return pd.DataFrame(keep)

######################################################################

def make_typemap(edges):
    nodes = pd.concat([edges["node1"], edges["node2"]])
    return pd.DataFrame({
        "n": nodes,
        "t": nodes.apply(lambda x: flags.get(x.split("__")[-1], "UNK"))
    }).drop_duplicates() #issues with deduplicaiting before, maybe this is the problem? 


def write_network(network, subnetworks):
    combined = pd.concat(subnetworks, ignore_index=True)
    combined.to_csv(f"intables/{network}_edges.csv", index=False)

    typemap = make_typemap(combined)
    typemap.to_csv(f"intables/{network}_typemap.csv", index=False, header=False)

# ====================================================================
pls_pls = pd.read_csv("edges/PLS-PLS.csv")
fec_fec = pd.read_csv("edges/FEC-FEC.csv")
pls_cpx = pd.read_csv("edges/PLS-CPX.csv")
fec_cpx = pd.read_csv("edges/FEC-CPX.csv")
fec_pls = pd.read_csv("edges/FEC-PLS.csv")

cc = pd.read_excel("edges/nodes 5.xlsx", sheet_name="Edges")

pls_pls_edges = load_edges(pls_pls, nodes)
print("PLS: ", len(pls_pls_edges))

fec_fec_edges = load_edges(fec_fec, nodes)
print("FEC: ", len(fec_fec_edges))

fec_pls_post = puc_fec_pls(fec_pls)
fec_pls_edges = load_edges(fec_pls_post, nodes, strip_suffix=True)
print("PLS-FEC post-puc: ", len(fec_pls_post))

pls_cpx_edges = load_edges(pls_cpx, nodes)
print("PLS-CPX: ", len(pls_cpx_edges))

fec_cpx_edges = load_edges(fec_cpx, nodes)
print("FEC-CPX: ", len(fec_cpx_edges))

cpx_cpx_edges = load_edges_excel(cc, "Choroid Plexus", "Choroid Plexus", nodes)
print("CPX-CPX: ", len(cpx_cpx_edges))
ctx_ctx_edges = load_edges_excel(cc, "Cortex", "Cortex", nodes, use_ids=True)
print("CTX-CTX:  ", len(ctx_ctx_edges))

stm_stm_edges = load_edges_excel(cc, "Striatum", "Striatum", nodes)
print("STM-STM: ", len(stm_stm_edges))

cpx_ctx_edges = load_edges_excel(cc, "Choroid Plexus", "Cortex", nodes, use_ids=True)
print("CPX-CTX: ", len(cpx_ctx_edges))

ctx_stm_edges = load_edges_excel(cc, "Cortex", "Striatum", nodes, use_ids=True)
print("CTX-STM: ", len(ctx_stm_edges))

networks = {
    "FEC-CTX": [
        fec_fec_edges,
        fec_pls_edges,
        pls_pls_edges,
        cpx_cpx_edges,
        pls_cpx_edges,
        cpx_ctx_edges,
        ctx_ctx_edges
    ],
    "FEC-STM_via_pls": [
        fec_fec_edges,
        fec_pls_edges,
        pls_pls_edges,
        cpx_cpx_edges,
        pls_cpx_edges,
        cpx_ctx_edges,
        ctx_stm_edges,
        stm_stm_edges
    ],
    "PLS-CTX": [
        pls_pls_edges,
        cpx_cpx_edges,
        cpx_ctx_edges,
        ctx_ctx_edges,
        pls_cpx_edges
    ],
    "PLS-STM": [
        pls_pls_edges,
        cpx_cpx_edges,
        cpx_ctx_edges,
        ctx_stm_edges,
        stm_stm_edges,
        pls_cpx_edges
    ],
    "FEC-CTX_via_cpx": [
        fec_fec_edges,
        fec_cpx_edges,
        cpx_cpx_edges,
        cpx_ctx_edges,
        ctx_ctx_edges
    ],
    "FEC-STM_via_cpx": [
        fec_fec_edges,
        fec_cpx_edges,
        cpx_cpx_edges,
        cpx_ctx_edges,
        ctx_stm_edges,
        stm_stm_edges
    ],
    "FEC-CTX_via_pls_cpx": [
        fec_fec_edges,
        fec_pls_edges,
        pls_pls_edges,
        fec_cpx_edges,
        cpx_cpx_edges,
        pls_cpx_edges,
        cpx_ctx_edges,
        ctx_ctx_edges
    ],
    "FEC-STM_via_pls_cpx": [
        fec_fec_edges,
        fec_pls_edges,
        pls_pls_edges,
        fec_cpx_edges,
        cpx_cpx_edges,
        pls_cpx_edges,
        cpx_ctx_edges,
        ctx_stm_edges,
        stm_stm_edges
    ],
    "PLS-CPX": [
        pls_pls_edges,
        cpx_cpx_edges,
        pls_cpx_edges
    ]
}


for network, subs in networks.items():
    write_network(network, subs)

print("done")
chime.success()
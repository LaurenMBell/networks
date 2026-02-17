import pandas as pd

def load_node(nodes:str, name:str, tissue, node_id):
    # Track every node’s tissue/ID so cross-tissue tables can infer labels.
    # nodes: {normalized_name: {"tissues": set(...), "ids": set(...)}}
    key = name.strip()
    entry = nodes.setdefault(key, {"tissues": set(), "ids": set()})
    entry["tissues"].add(tissue)
    if node_id is not None:
        entry["ids"].add(str(node_id))

def node_tissue(nodes, name: str):
    # Returns the tissue if we’ve only seen this node in one context.
    entry = nodes.get(name.strip())
    if not entry:
        return None
    tissues = entry["tissues"]
    return next(iter(tissues)) if len(tissues) == 1 else None

def canonical_name(name: str, tissue: str, node_id=None) -> str:
    # Normalize names and tag CTX “UNKNOWN” nodes with their IDs so they’re unique.
    cleaned = str(name).strip()
    if cleaned == "UNKNOWN" and tissue == "Cortex" and node_id is not None:
        return f"{cleaned}_{node_id}"
    return cleaned


def label_node(name:str, tissue):
    # Append the __CPX/PLS/etc suffix that other scripts expect.
    flag = tissue_flag[tissue]
    suffix = f"__{flag}"
    return name if name.endswith(suffix) else f"{name}{suffix}"


def strip_flag(text: str) -> str:
    # Remove “-P/-F” suffixes that appear in the raw FEC-PLS tables.
    if text.endswith("-P") or text.endswith("-F"):
        return text[:-2]
    return text


def load_edges(df: pd.DataFrame, nodes) -> None:
    # Feed the single-tissue edge lists so later tables can inherit tissues.
    for i, row in df.iterrows():
        load_node(nodes, canonical_name(row["Node 1 Name"], row["Node 1 Tissue"]),
            row["Node 1 Tissue"])
        load_node(nodes, canonical_name(row["Node 2 Name"], row["Node 2 Tissue"]),
            row["Node 2 Tissue"])


def load_cpx_ctx_edges(cc: pd.DataFrame, nodes) -> pd.DataFrame:
    allowed = {"Choroid Plexus", "Cortex"}
    mask = cc["Node 1 Tissue"].isin(allowed) & cc["Node 2 Tissue"].isin(allowed)
    subset = cc.loc[mask]
    edges = []
    for i, row in subset.iterrows():
        t1 = row["Node 1 Tissue"].strip()
        t2 = row["Node 2 Tissue"].strip()
        if t1 == t2:
            continue
        name1 = canonical_name(row["Node 1 Name"], t1, row.get("Node 1 ID"))
        name2 = canonical_name(row["Node 2 Name"], t2, row.get("Node 2 ID"))
        load_node(nodes, name1, t1, row.get("Node 1 ID"))
        load_node(nodes, name2, t2, row.get("Node 2 ID"))
        edges.append({
            "node1": label_node(name1, t1),
            "node2": label_node(name2, t2),
        })
    return pd.DataFrame(edges)


def load_pls_cpx_edges(pls_cpx: pd.DataFrame, nodes) -> pd.DataFrame:
    edges = []
    for i, row in pls_cpx.iterrows():
        raw_cpx = str(row["n1"]).strip()
        tissue_cpx = node_tissue(nodes, raw_cpx)
        load_node(nodes, raw_cpx, tissue_cpx)

        raw_pls = str(row["n2"]).strip()
        tissue_pls = node_tissue(nodes, raw_pls)
        load_node(nodes, raw_pls, tissue_pls)

        edges.append({"node1": label_node(raw_cpx, tissue_cpx),
                    "node2": label_node(raw_pls, tissue_pls)})
    return pd.DataFrame(edges)


def puc_fec_pls(df, threshold) -> pd.DataFrame:
    
    adjacency = {}
    dirs = {}
    for idx, row in df.iterrows():
        direction = int(row["edge_dir"])
        dirs[idx] = direction
        adjacency.setdefault(str(row["Node 1 Name"]).strip(), []).append(idx)
        adjacency.setdefault(str(row["Node 2 Name"]).strip(), []).append(idx)

    to_remove = set()
    for indices in adjacency.values():
        pos = sum(1 for idx in indices if dirs[idx] > 0)
        neg = sum(1 for idx in indices if dirs[idx] < 0)
        total = pos + neg
        if total == 0:
            continue
        minority = min(pos, neg)
        if minority == 0:
            continue
        ratio = minority / total
        if ratio <= threshold:
            majority_sign = 1 if pos >= neg else -1
            for idx in indices:
                if dirs[idx] != majority_sign:
                    to_remove.add(idx)
        else:
            to_remove.update(indices)

    if not to_remove:
        return df.reset_index(drop=True)
    return df.drop(index=to_remove).reset_index(drop=True)


def load_fec_pls_edges(fec_pls, nodes):
    edges = []
    for i, row in fec_pls.iterrows():
        fec_name = canonical_name(strip_flag(row["Node 1 Name"]), "Feces")
        load_node(nodes, fec_name, "Feces")

        tissue_pls = row.get("Node 2 Tissue", "Plasma")
        tissue_pls = str(tissue_pls).strip() if not pd.isna(tissue_pls) else "Plasma"
        if tissue_pls not in tissue_flag:
            tissue_pls = "Plasma"
        pls_name = canonical_name(strip_flag(row["Node 2 Name"]), tissue_pls)
        load_node(nodes, pls_name, tissue_pls)

        edges.append({
            "node1": label_node(fec_name, "Feces"),
            "node2": label_node(pls_name, tissue_pls),
        })
    return pd.DataFrame(edges)


def extract_tissue(node:str):
    return flag_map.get(node.split("__")[-1], "UNK")


def make_typemap(edges):
    nodes = pd.concat([edges["node1"], edges["node2"]]).drop_duplicates().reset_index(drop=True)
    return pd.DataFrame({"n": nodes, "t": nodes.apply(extract_tissue)})


def write_network(network, subnetworks):
    combined = pd.concat(subnetworks, ignore_index=True)

    combined.to_csv(f"intables/{network}_edges.csv", index=False)

    typemap = make_typemap(combined)
    typemap.to_csv(f"intables/{network}_typemap.csv", index=False, header=False)


#====================================================================================================

#DICT WITH ALL NODES : tissues/id
#{node : {tissues: set(), ids: set()}}
nodes = {}

tissue_flag = {"Plasma": "PLS", "Feces": "FEC",
    "Choroid Plexus": "CPX","Cortex": "CTX"}

flag_map = {v: k for k, v in tissue_flag.items()}

pls_truth = pd.read_csv("edges/PLS-PLS.csv")
load_edges(pls_truth, nodes)

fec_truth = pd.read_csv("edges/FEC-FEC.csv")
load_edges(fec_truth, nodes)

cc = pd.read_excel("edges/nodes 4.xlsx", sheet_name="Edges")
cpx_ctx_edges = load_cpx_ctx_edges(cc, nodes)

pls_cpx_df = pd.read_csv("edges/PLS-CPX.csv")
pls_cpx_edges = load_pls_cpx_edges(pls_cpx_df, nodes)

fec_pls = pd.read_csv("edges/FEC-PLS_prerpuc.csv")
fec_pls_cleaned = puc_fec_pls(fec_pls, 0.2)
fec_pls_edges = load_fec_pls_edges(fec_pls_cleaned, nodes)

networks = {"PLS-CTX": [cpx_ctx_edges, pls_cpx_edges],
    "FEC-CTX": [cpx_ctx_edges, pls_cpx_edges, fec_pls_edges]}

for network, subs in networks.items():
    write_network(network, subs)

print("done with typemap")

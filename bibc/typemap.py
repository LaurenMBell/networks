import pandas as pd
import pickle
import numpy as np

def load_node(nodes, name, tissue, node_id=None):
    #{name: {tissues: set()), ids: set()}}
    entry = nodes.setdefault(name.strip(), {"tissues": set(), "ids": set()})
    entry["tissues"].add(tissue)
    if node_id is not None:
        entry["ids"].add(str(node_id))


def name_nodes(name, tissue, node_id=None):
    #strip node names and find CTX UNKNOWNS
    cleaned = str(name).strip()
    if (
        cleaned == "UNKNOWN"
        and tissue in {"Cortex", "Striatum"}
        and node_id is not None
    ):
        return f"{cleaned}_{node_id}"
    return cleaned


def label_node(name, tissue):
    flag = tissue_flag[tissue]
    tag = f"__{flag}"
    if name.endswith(tag):
        return name
    return f"{name}{tag}"


def load_cpx_ctx_edges(cc, nodes):
    df = cc.rename(columns={"Node 1 Name": "node1", "Node 2 Name": "node2"}).copy()
    df["node1"] = df["node1"].apply(lambda v: "UNKNOWN" if str(v).strip() == "UNKNOWN" else v)
    df["node2"] = np.where(
        df["node2"].eq("UNKNOWN"),
        "UNKNOWN_" + df["Node 2 ID"].astype(str),
        df["node2"],
    )

    mask = (df["Node 1 Tissue"] == "Choroid Plexus") & (df["Node 2 Tissue"] == "Cortex")
    df = df.loc[mask, [
        "node1",
        "node2",
        "Node 1 ID",
        "Node 2 ID",
        "Node 1 Tissue",
        "Node 2 Tissue",
    ]].copy()

    if df.empty:
        return pd.DataFrame(columns=["node1", "node2"])

    df["node_min"] = df[["node1", "node2"]].min(axis=1)
    df["node_max"] = df[["node1", "node2"]].max(axis=1)
    df = df.drop_duplicates(subset=["node_min", "node_max"])
    df = df.drop(columns=["node_min", "node_max"])
    df = df.apply(lambda col: col.str.strip() if col.dtype == "object" else col)

    edges = []
    for _, row in df.iterrows():
        name1 = row["node1"]
        name2 = row["node2"]

        load_node(nodes, name1, "Choroid Plexus", row.get("Node 1 ID"))
        load_node(nodes, name2, "Cortex", row.get("Node 2 ID"))

        edges.append({
            "node1": label_node(name1, "Choroid Plexus"),
            "node2": label_node(name2, "Cortex"),
        })

    return pd.DataFrame(edges)

def load_ctx_stm_edges(cc, nodes):
    df = cc.rename(columns={"Node 1 Name": "node1", "Node 2 Name": "node2"}).copy()
    mask = (
        (df["Node 1 Tissue"] == "Cortex") & (df["Node 2 Tissue"] == "Striatum")
    ) | (
        (df["Node 1 Tissue"] == "Striatum") & (df["Node 2 Tissue"] == "Cortex")
    )
    df = df.loc[mask, [
        "node1",
        "node2",
        "Node 1 ID",
        "Node 2 ID",
        "Node 1 Tissue",
        "Node 2 Tissue",
    ]].copy()

    if df.empty:
        return pd.DataFrame(columns=["node1", "node2"])

    df = df.apply(lambda col: col.str.strip() if col.dtype == "object" else col)

    df["ctx_name"] = np.where(df["Node 1 Tissue"] == "Cortex", df["node1"], df["node2"])
    df["ctx_id"] = np.where(df["Node 1 Tissue"] == "Cortex", df["Node 1 ID"], df["Node 2 ID"])
    df["stm_name"] = np.where(df["Node 1 Tissue"] == "Striatum", df["node1"], df["node2"])
    df["stm_id"] = np.where(df["Node 1 Tissue"] == "Striatum", df["Node 1 ID"], df["Node 2 ID"])

    df["ctx_name"] = [
        name_nodes(name, "Cortex", node_id)
        for name, node_id in zip(df["ctx_name"], df["ctx_id"])
    ]
    df["stm_name"] = [
        name_nodes(name, "Striatum", node_id)
        for name, node_id in zip(df["stm_name"], df["stm_id"])
    ]

    df = df.drop_duplicates(subset=["ctx_name", "stm_name"])

    edges = []
    for _, row in df.iterrows():
        ctx_name = row["ctx_name"]
        stm_name = row["stm_name"]
        ctx_id = row["ctx_id"]
        stm_id = row["stm_id"]

        load_node(nodes, ctx_name, "Cortex", ctx_id)
        load_node(nodes, stm_name, "Striatum", stm_id)

        edges.append({
            "node1": label_node(ctx_name, "Cortex"),
            "node2": label_node(stm_name, "Striatum"),
        })

    print("ctx-str")

    return pd.DataFrame(edges)

def dedupe_undirected(df):
    if df.empty:
        return df
    df = df.copy()
    df["node_min"] = df[["node1", "node2"]].min(axis=1)
    df["node_max"] = df[["node1", "node2"]].max(axis=1)
    df = df.drop_duplicates(subset=["node_min", "node_max"])
    return df.drop(columns=["node_min", "node_max"])


def load_cpx_cpx_edges(cc, nodes):
    df = cc.rename(columns={"Node 1 Name": "node1", "Node 2 Name": "node2"}).copy()
    df["node1"] = df["node1"].apply(lambda v: "UNKNOWN" if str(v).strip() == "UNKNOWN" else v)
    df["node2"] = df["node2"].apply(lambda v: "UNKNOWN" if str(v).strip() == "UNKNOWN" else v)

    mask = (df["Node 1 Tissue"] == "Choroid Plexus") & (df["Node 2 Tissue"] == "Choroid Plexus")
    df = df.loc[mask, ["node1", "node2", "Node 1 ID", "Node 2 ID"]]
    df = df.apply(lambda col: col.str.strip() if col.dtype == "object" else col)
    df = dedupe_undirected(df)

    edges = []
    for _, row in df.iterrows():
        name1 = row["node1"]
        name2 = row["node2"]

        load_node(nodes, name1, "Choroid Plexus", row.get("Node 1 ID"))
        load_node(nodes, name2, "Choroid Plexus", row.get("Node 2 ID"))

        edges.append({
            "node1": label_node(name1, "Choroid Plexus"),
            "node2": label_node(name2, "Choroid Plexus"),
        })

    return pd.DataFrame(edges)


def load_ctx_ctx_edges(cc, nodes):
    df = cc.rename(columns={"Node 1 Name": "node1", "Node 2 Name": "node2"}).copy()
    df["node1"] = np.where(
        df["node1"].eq("UNKNOWN"),
        "UNKNOWN_" + df["Node 1 ID"].astype(str),
        df["node1"],
    )
    df["node2"] = np.where(
        df["node2"].eq("UNKNOWN"),
        "UNKNOWN_" + df["Node 2 ID"].astype(str),
        df["node2"],
    )

    mask = (df["Node 1 Tissue"] == "Cortex") & (df["Node 2 Tissue"] == "Cortex")
    df = df.loc[mask, ["node1", "node2", "Node 1 ID", "Node 2 ID"]]
    df = df.apply(lambda col: col.str.strip() if col.dtype == "object" else col)
    df = dedupe_undirected(df)

    edges = []
    for i, row in df.iterrows():
        name1 = row["node1"]
        name2 = row["node2"]

        load_node(nodes, name1, "Cortex", row.get("Node 1 ID"))
        load_node(nodes, name2, "Cortex", row.get("Node 2 ID"))

        edges.append({
            "node1": label_node(name1, "Cortex"),
            "node2": label_node(name2, "Cortex"),
        })

    return pd.DataFrame(edges)

def load_stm_stm_edges(cc, nodes):
    df = cc.rename(columns={"Node 1 Name": "node1", "Node 2 Name": "node2"})

    mask = (df["Node 1 Tissue"] == "Striatum") & (df["Node 2 Tissue"] == "Striatum")
    df = df.loc[mask, ["node1", "node2", "Node 1 ID", "Node 2 ID"]]
    df = df.apply(lambda col: col.str.strip() if col.dtype == "object" else col)
    df = dedupe_undirected(df)

    edges = []
    for i, row in df.iterrows():
        name1 = row["node1"]
        name2 = row["node2"]

        load_node(nodes, name1, "Striatum", row.get("Node 1 ID"))
        load_node(nodes, name2, "Striatum", row.get("Node 2 ID"))

        edges.append({
            "node1": label_node(name1, "Striatum"),
            "node2": label_node(name2, "Striatum")})

    return pd.DataFrame(edges)


def load_pls_cpx_edges(pls_cpx, nodes):
    edges = []
    for i, row in pls_cpx.iterrows():
        t1 = row["Node 1 Tissue"]
        t2 = row["Node 2 Tissue"]

        name1 = name_nodes(row["Node 1 Name"], t1)
        name2 = name_nodes(row["Node 2 Name"], t2)

        load_node(nodes, name1, t1)
        load_node(nodes, name2, t2)

        edges.append({"node1": label_node(name1, t1),
                    "node2": label_node(name2, t2)})
    return pd.DataFrame(edges)

def load_fec_cpx_edges(fec_cpx, nodes):
    edges = []
    for i, row in fec_cpx.iterrows():
        t1 = row["Node 1 Tissue"]
        t2 = row["Node 2 Tissue"]

        name1 = name_nodes(row["Node 1 Name"], t1)
        name2 = name_nodes(row["Node 2 Name"], t2)

        load_node(nodes, name1, t1)
        load_node(nodes, name2, t2)

        edges.append({"node1": label_node(name1, t1),
                    "node2": label_node(name2, t2)})
    return pd.DataFrame(edges)

def strip_tissue_suffix(raw_name):
    value = str(raw_name).strip()
    if value.endswith("-F") or value.endswith("-P"):
        return value[:-2]
    return value


def build_direction_map(df):
    mapping = {}
    for i, row in df.iterrows():
        name = strip_tissue_suffix(row["node"])
        mapping[name] = row["node_dir"]
    return mapping


def puc_fec_pls(df):
    fec_df = pd.read_csv("edges/FEC_nodes.csv")
    pls_df = pd.read_csv("edges/PLS_nodes.csv")

    fec_node_dir = build_direction_map(fec_df)
    pls_node_dir = build_direction_map(pls_df)

    to_remove = []

    for idx, row in df.iterrows():
        pls_name = strip_tissue_suffix(row["Node 1 Name"])
        fec_name = strip_tissue_suffix(row["Node 2 Name"])
        edge_dir = row["edge_dir"]

        dir_pls = pls_node_dir.get(pls_name)
        dir_fec = fec_node_dir.get(fec_name)

        if edge_dir is None or dir_pls is None or dir_fec is None:
            #print(f"NO DIR: PLS {pls_name} ({dir_pls}) or FEC {fec_name} ({dir_fec})")
            to_remove.append(idx)
            continue

        if edge_dir * dir_fec != dir_pls:
            to_remove.append(idx)

    cleaned = df.drop(index=to_remove).reset_index(drop=True)
    return cleaned


def load_fec_pls_edges(fec_pls, nodes):
    edges = []
    for i, row in fec_pls.iterrows():
        pls_name = name_nodes(strip_tissue_suffix(row["Node 1 Name"]), "Plasma")
        
        load_node(nodes, pls_name, "Plasma")

        fec_name = name_nodes(strip_tissue_suffix(row["Node 2 Name"]), "Feces")
        load_node(nodes, fec_name, "Feces")

        edges.append({
            "node1": label_node(fec_name, "Feces"),
            "node2": label_node(pls_name, "Plasma")
        })
    return pd.DataFrame(edges)


def load_edges(df, tissue, nodes):
    edges = []
    for i, row in df.iterrows():
        name1 = name_nodes(row["Node 1 Name"], tissue, row.get("Node 1 ID"))
        name2 = name_nodes(row["Node 2 Name"], tissue, row.get("Node 2 ID"))

        load_node(nodes, name1, tissue, row.get("Node 1 ID"))
        load_node(nodes, name2, tissue, row.get("Node 2 ID"))

        edges.append({"node1": label_node(name1, tissue),
                      "node2": label_node(name2, tissue)})
    return pd.DataFrame(edges) if edges else pd.DataFrame(columns=["node1", "node2"])


def make_typemap(edges):
    nodes = pd.concat([edges["node1"], edges["node2"]]).reset_index(drop=True)
    return pd.DataFrame({
        "n": nodes,
        "t": nodes.apply(lambda x: flag_map.get(x.split("__")[-1], "UNK"))})


def write_network(network, subnetworks):
    combined = pd.concat(subnetworks, ignore_index=True)

    combined.to_csv(f"intables/{network}_edges.csv", index=False)

    typemap = make_typemap(combined)
    typemap.drop_duplicates().to_csv(f"intables/{network}_typemap.csv", index=False, header=False)


#====================================================================================================

#DICT WITH ALL NODES : tissues/id
#{node : {tissues: set(), ids: set()}}
nodes = {}

tissue_flag = {
    "Plasma": "PLS",
    "Feces": "FEC",
    "Choroid Plexus": "CPX",
    "Cortex": "CTX",
    "Striatum": "STM",
}

flag_map = {v: k for k, v in tissue_flag.items()}

pls_truth = pd.read_csv("edges/PLS-PLS.csv")
pls_pls_edges = load_edges(pls_truth, "Plasma", nodes)
print("PLS: ", len(pls_pls_edges))

fec_truth = pd.read_csv("edges/FEC-FEC.csv")
fec_fec_edges = load_edges(fec_truth, "Feces", nodes)
print("FEC: ", len(fec_fec_edges))

fec_cpx = pd.read_csv("edges/FEC-CPX.csv")
fec_cpx_edges = load_fec_cpx_edges(fec_cpx, nodes)
print("FEC-CPX: ", len(fec_cpx_edges))

cc = pd.read_excel("edges/nodes 5.xlsx", sheet_name="Edges")
cpx_ctx_edges = load_cpx_ctx_edges(cc, nodes)
print("CPX-CTX: ", len(cpx_ctx_edges))
cpx_cpx_edges = load_cpx_cpx_edges(cc, nodes)
print("CPX-CPX: ", len(cpx_cpx_edges))
ctx_ctx_edges = load_ctx_ctx_edges(cc, nodes)
print("CTX-CTX:  ", len(ctx_ctx_edges))
ctx_stm_edges = load_ctx_stm_edges(cc, nodes)
print("CTX-STM: ", len(ctx_stm_edges))
stm_stm_edges = load_stm_stm_edges(cc, nodes)
print("STM-STM: ", len(stm_stm_edges))

pls_cpx_df = pd.read_csv("edges/PLS-CPX.csv")
pls_cpx_edges = load_pls_cpx_edges(pls_cpx_df, nodes)
print("PLS-CPX: ", len(pls_cpx_edges))

fec_pls = pd.read_csv("edges/FEC-PLS.csv")
fec_pls_post = puc_fec_pls(fec_pls)
print("PLS-FEC post-puc: ", len(fec_pls_post))
fec_pls_edges = load_fec_pls_edges(fec_pls_post, nodes)

networks = {
    "FEC-CTX": [
        fec_fec_edges,
        fec_pls_edges,
        pls_pls_edges,
        cpx_cpx_edges,
        pls_cpx_edges,
        cpx_ctx_edges,
        ctx_ctx_edges,
    ],
    "FEC-STM": [
        fec_fec_edges,
        fec_pls_edges,
        pls_pls_edges,
        cpx_cpx_edges,
        pls_cpx_edges,
        cpx_ctx_edges,
        ctx_stm_edges,
        stm_stm_edges,
    ],
    "PLS-CTX": [
        pls_pls_edges,
        cpx_cpx_edges,
        cpx_ctx_edges,
        ctx_ctx_edges,
        pls_cpx_edges,
    ],
    "PLS-STM": [
        pls_pls_edges,
        cpx_cpx_edges,
        cpx_ctx_edges,
        ctx_stm_edges,
        stm_stm_edges,
        pls_cpx_edges,
    ],
    "FEC-CTX_via_cpx": [
        fec_fec_edges,
        fec_cpx_edges,
        cpx_cpx_edges,
        cpx_ctx_edges,
        ctx_ctx_edges,
    ],
    "FEC-STM_via_cpx": [
        fec_fec_edges,
        fec_cpx_edges,
        cpx_cpx_edges,
        cpx_ctx_edges,
        ctx_stm_edges,
        stm_stm_edges,
    ],
    "FEC-CTX_via_pls_cpx": [
        fec_fec_edges,
        fec_pls_edges,
        pls_pls_edges,
        fec_cpx_edges,
        cpx_cpx_edges,
        pls_cpx_edges,
        cpx_ctx_edges,
        ctx_ctx_edges,
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
        stm_stm_edges,
    ],
    "PLS-CPX": [
        pls_pls_edges,
        cpx_cpx_edges,
        pls_cpx_edges,
    ],
}

for network, subs in networks.items():
    write_network(network, subs)

print("done with typemap")

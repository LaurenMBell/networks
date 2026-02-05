import pandas as pd

#deal with UNKNOWN genes in ctx
def unknowns(name, tissue, node_id):
    if name == "UNKNOWN" and tissue == "Cortex":
        return f"{name}_{node_id}"
    else:
        return name
    
#make typemap - HOW TO ISOLATE UNIQUE EDGES 
def filter_edges(df : pd.DataFrame, valid_pairs):
    edges = df[df.apply(lambda r: (r["Node 1 Tissue"], r["Node 2 Tissue"]) in valid_pairs, axis=1)]

    """
    edge_lst = set()
    for edge in edges.rows:
        e =  [edge["Node 1"], edge["Node 2"]]
        edge_lst.add(e) """

    

    return edges

def extract_tissue(node):

    flag = node.split("__")[-1]
    return flag_map.get(flag, "UNK")

def build_typemap(edge_file, out_file):
    edges = pd.read_csv(edge_file)
    nodes = pd.concat([edges["node1"], edges["node2"]]).drop_duplicates().reset_index(drop=True)

    type_map = pd.DataFrame({
        "n" : nodes,
        "t" :nodes.apply(extract_tissue)})

    type_map.to_csv(out_file, index=False, header=False)
    


df = pd.read_csv("matt_edges.csv")

tissue_flag = {"Plasma": "PLS", "Feces": "FEC",
    "Choroid Plexus": "CPX","Cortex": "CTX"}

flag_map = {v: k for k, v in tissue_flag.items()}

df = df[~df["Node 1 Tissue"].isin(["Remaining Brain", "Striatum"])]
df = df[~df["Node 2 Tissue"].isin(["Remaining Brain", "Striatum"])]

df["Node1_unique"] = df.apply(
    lambda r: unknowns(r["Node 1 Name"], r["Node 1 Tissue"], r["Node 1 ID"]), axis=1)
df["Node2_unique"] = df.apply(
    lambda r: unknowns(r["Node 2 Name"], r["Node 2 Tissue"], r["Node 2 ID"]), axis=1)
df["Node1_labeled"] = df["Node1_unique"] + "__" + df["Node 1 Tissue"].map(tissue_flag)
df["Node2_labeled"] = df["Node2_unique"] + "__" + df["Node 2 Tissue"].map(tissue_flag)


pls_ctx_pairs = {("Plasma", "Plasma"), ("Plasma", "Choroid Plexus"), ("Choroid Plexus", "Choroid Plexus"),
                  ("Choroid Plexus", "Plasma"), ("Choroid Plexus", "Cortex"), ("Cortex", "Cortex"), ("Cortex", "Choroid Plexus")}

fec_ctx_pairs = {("Feces","Feces"), ("Feces", "Plasma"),("Plasma", "Plasma"), ("Plasma", "Feces"),("Plasma", "Choroid Plexus"),
    ("Choroid Plexus", "Choroid Plexus"), ("Choroid Plexus", "Plasma"),("Choroid Plexus", "Cortex"),("Cortex", "Cortex"),
    ("Cortex", "Choroid Plexus")}

#pls-ctx
pls_ctx = filter_edges(df, pls_ctx_pairs)

pls_ctx_edges = pls_ctx[["Node1_labeled", "Node2_labeled"]]
pls_ctx_edges.columns = ["node1", "node2"]
pls_ctx_edges.to_csv("PLS-CTX_edges.csv", index=False)

build_typemap("PLS-CTX_edges.csv", "PLS-CTX_typemap.csv")

#fec-ctx
fec_ctx = filter_edges(df, fec_ctx_pairs)

fec_ctx_edges = fec_ctx[["Node1_labeled", "Node2_labeled"]]
fec_ctx_edges.columns = ["node1", "node2"]
fec_ctx_edges.to_csv("FEC-CTX_edges.csv", index=False)

build_typemap("FEC-CTX_edges.csv", "FEC-CTX_typemap.csv")

print("done with typemap")

## STEP 1.7 OF ANALYSIS PLAN

import pandas as pd
import pickle
import numpy as np

p = open("../id_to_symbol_map.pickle", 'rb')
name_dict = pickle.load(p)
p.close()

def predict_one_sign(row):
        csf = row["csf_dir"]
        brain = row["brain_dir"]
        if pd.isna(csf):
            return np.nan
        if pd.isna(brain):
             return np.nan
        if csf == brain:
            return "1"
        else:
            return "-1"
        
def predict_signs(edges: pd.DataFrame, brain_nodes: pd.DataFrame):
    df = edges.copy()
    
    brain_dirs = brain_nodes.set_index("ID")["Mean Log2 Fold Change Direction (DSS)"]
    df["brain_dir"] = df["n1"].map(brain_dirs)

    csf_dirs = csf_nodes.set_index("node")["node_dir"]
    df["csf_dir"] = df["n2"].map(csf_dirs)

    df["predicted_sign"] = df.apply(predict_one_sign, axis=1)
    return df[["n1", "n2", "csf_dir", "brain_dir", "predicted_sign"]]

#p = open("../id_to_symbol_map.pickle", 'rb')
#name_dict = pickle.load(p)
#p.close()

csf_nodes = pd.read_csv("../csf_rpuc_nodes.csv")
brain_node_dir = pd.read_csv("../network_nodes.csv")

ctx_node_dir = brain_node_dir[brain_node_dir["Type"] == "gene_C"]
rbn_node_dir = brain_node_dir[brain_node_dir["Type"] == "gene_B"]
str_node_dir = brain_node_dir[brain_node_dir["Type"] == "gene_S"]

csf_ctx_edges = pd.read_csv("../../CSF-CTX/CSF-CTX_FDR_threshold.csv")
csf_rbn_edges = pd.read_csv("../../CSF-RBN/CSF-RBN_FDR_threshold.csv")
csf_str_edges = pd.read_csv("../../CSF-STR/CSF-STR_FDR_threshold.csv")

ctx_predicted = predict_signs(csf_ctx_edges, ctx_node_dir)
rbn_predicted = predict_signs(csf_rbn_edges, rbn_node_dir)
str_predicted = predict_signs(csf_str_edges, str_node_dir)

ctx_predicted.to_csv("csf_ctx_predicted_signs.csv", index=False)
rbn_predicted.to_csv("csf_rbn_predicted_signs.csv", index=False)
str_predicted.to_csv("csf_str_predicted_signs.csv", index=False)
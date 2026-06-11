"""
import pandas as pd 

csf_str_edges = pd.read_csv("../CSF-STR/CSF-STR_FDR_threshold.csv")
csf_rbn_edges = pd.read_csv("../CSF-RBN/CSF-RBN_FDR_threshold.csv")
csf_ctx_edges = pd.read_csv("../CSF-CTX/CSF-CTX_FDR_threshold.csv")

csf_str_pred_signs = pd.read_csv("../CSF_ESP/CSF-brain_edge_signs/csf_str_predicted_signs.csv")
csf_rbn_pred_signs = pd.read_csv("../CSF_ESP/CSF-brain_edge_signs/csf_rbn_predicted_signs.csv")
csf_ctx_pred_signs = pd.read_csv("../CSF_ESP/CSF-brain_edge_signs/csf_ctx_predicted_signs.csv")

final_table = pd.DataFrame([["n1", "n2", "VECPAC_r", "LPS_r", "DSS_r", "pooled_r", 
             "pred_sign", "model_agreement_pooled_sign", 
             "model_agreement_pred_sign", "pooled_pred_agreement_sign", 
             "all_sign_agreement"]])
"""

import pandas as pd

csf_str_edges = pd.read_csv("../CSF-STR/CSF-STR_FDR_threshold.csv")
csf_rbn_edges = pd.read_csv("../CSF-RBN/CSF-RBN_FDR_threshold.csv")
csf_ctx_edges = pd.read_csv("../CSF-CTX/CSF-CTX_FDR_threshold.csv")
csf_str_pred_signs = pd.read_csv("../CSF_ESP/CSF-brain_edge_signs/csf_str_predicted_signs.csv")
csf_rbn_pred_signs = pd.read_csv("../CSF_ESP/CSF-brain_edge_signs/csf_rbn_predicted_signs.csv")
csf_ctx_pred_signs = pd.read_csv("../CSF_ESP/CSF-brain_edge_signs/csf_ctx_predicted_signs.csv")

csf_str_edges["subnetwork"] = "STR"
csf_rbn_edges["subnetwork"] = "RBN"
csf_ctx_edges["subnetwork"] = "CTX"
csf_str_pred_signs["subnetwork"] = "STR"
csf_rbn_pred_signs["subnetwork"] = "RBN"
csf_ctx_pred_signs["subnetwork"] = "CTX"

edges = pd.concat([csf_str_edges, csf_rbn_edges, csf_ctx_edges])
pred_signs = pd.concat([csf_str_pred_signs, csf_rbn_pred_signs, csf_ctx_pred_signs]) 

#lauren got help from AI for this step
df = edges[["n1", "n2", "subnetwork", "VECPAC_r", "LPS_r", "DSS_r", "Pooled_r", "Sign"]].merge(
    pred_signs[["n1", "n2", "subnetwork", "predicted_sign"]],
    on=["n1", "n2", "subnetwork"],
    how="left")

df["pooled_r"] = df["Pooled_r"]
df["model_sign_num"] = df["Sign"].map({"+": 1, "-": -1})
df["pooled_r_sign"] = df["pooled_r"].apply(lambda x: 1 if x > 0 else -1)



df["model_agreement_pooled_sign"] = (df["model_sign_num"] == df["pooled_r_sign"]).map({True: 1, False: -1})
df["model_agreement_pred_sign"] = (df["model_sign_num"] == df["predicted_sign"]).map({True: 1, False: -1})
df["pooled_pred_agreement_sign"] = (df["pooled_r_sign"] == df["predicted_sign"]).map({True: 1, False: -1})

df["all_sign_agreement"] = ((df["model_agreement_pooled_sign"] == 1) &
    (df["model_agreement_pred_sign"] == 1) & 
    (df["pooled_pred_agreement_sign"] == 1)).map({True: 1, False: -1})

final_table = df[["n1", "n2", "VECPAC_r", "LPS_r", "DSS_r", "pooled_r",
                   "predicted_sign", "model_agreement_pooled_sign",
                   "model_agreement_pred_sign", "pooled_pred_agreement_sign",
                   "all_sign_agreement"]]

final_table.to_csv("CSF_analysis_table.csv", index=False)

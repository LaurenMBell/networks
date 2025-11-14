import pandas as pd
import numpy as np

# FUNCS TO CALC NETWORK CHARACTERISTICS FOR EACH PARAM ####################

def filter_thresholds(ind, fdr, df):
    #takes total fdr table, returns filtered one

    filtered_df = df[df["FDR"] <= fdr]
    max_p_val = filtered_df[["VECPAC p-values", "DSS p-values", "LPS p-values"]].max(axis=1)
    filtered_df = filtered_df[max_p_val <= ind]

    return filtered_df #new dataframe

def nodes(df):
    # takes filtered df, counts uniqe metabolites, then finds + and - mets
    
    #to be changed later
    pos_nodes = np.nan
    neg_nodes = np.nan

    unique_mets = set() 

    unique_mets = set(df["Metabolite 1"]).union(set(df["Metabolite 2"]))
    total_nodes = len(unique_mets)
    
    # question: how to calculate positive and negative nodes, fix this lauren
    return total_nodes, pos_nodes, neg_nodes 

def edges(df):
    # takes filtered df, counts number of corr coefs
    #finds + and - correlations and the ratio btwn the two

    total_edges = len(df)

    to_check = ["VECPAC r", "LPS r", "DSS r"]

    #nonnan count is a row with true or false if it is or isn't na
    #so adding them up should be 2 or 3
    non_nan_count = df[to_check].notna().sum(axis=1)  # Series
    sign_sum = np.nansum(np.sign(df[to_check].values), axis=1)  # Series

    # Element-wise comparisons create boolean masks
    pos_edges = (sign_sum == non_nan_count).sum()  # Count True values
    neg_edges = (sign_sum == -non_nan_count).sum()

    #pos_edges = len(df[df["edge_dir"] == 1])
    #neg_edges = len(df[df["edge_dir"] == -1])

    if neg_edges > 0:
        ratio = pos_edges/neg_edges 
    else:
        ratio = np.nan

    return total_edges, pos_edges, neg_edges, ratio

def mean_corr(df):
    #take mean of each edge, then mean of those measn 
    
    # Only take correlation columns
    r_cols = ["VECPAC r", "DSS r", "LPS r"]

    edge_means = df[r_cols].apply(pd.to_numeric, errors="coerce").mean(axis=1)
    overall_mean = edge_means.mean()

    return overall_mean


# constructs the final table for each file !!
def table_making(file, name):
    final_table = pd.DataFrame()

    #corr_cols = ["VECPAC r", "LPS r", "DSS r"]
    #full_graph = file.dropna(subset=corr_cols)  
    
    full_graph = file.copy()  
    total_edges_full = len(full_graph)

    for i in ind_p:
        for j in fdr_p:
            return_table = pd.DataFrame()
            
            for col in ["VECPAC r", "VECPAC p-values", "DSS r", "DSS p-values", "LPS r", "LPS p-values", "FDR"]:
                file[col] = pd.to_numeric(file[col], errors="coerce")
            file.columns = file.columns.str.strip()

        
            filtered = filter_thresholds(i, j, full_graph)
            #print(len(filtered))


            total_nodes, pos_nodes, neg_nodes = nodes(filtered)
            total_edges, pos_edges, neg_edges, ratio = edges(filtered)
            full_edges = (total_nodes**2 - total_nodes)/2
            density = total_edges/full_edges
            mean_corr_coeff = mean_corr(filtered)

            return_table = pd.DataFrame([{
                "Individual_Pval": i,
                "FDR": j,
                "Total_Nodes": total_nodes,
                #"Positive_Nodes": pos_nodes,
                #"Negative_Nodes": neg_nodes,
                "Total_Edges": total_edges,
                "Positive_Edges": pos_edges,
                "Negative_Edges": neg_edges,
                "Pos_to_Neg": ratio,
                "Density": density,
                "Mean_Corr_Coeff": mean_corr_coeff
            }])

            final_table = pd.concat([final_table, return_table])

    final_table.to_csv(f"{name}_network_properties.csv", index=False)

    return final_table


# TABLE MAKING #########################################################################

#df_quant = pd.read_csv("/Users/laurenbell/Desktop/Metabolite_Brain_Gut_MorgunLab/network_properties/quantile_norm_data/updated_metaanalysis_fdr_table.csv")
#df_median_with_zeros = pd.read_csv("/Users/laurenbell/Desktop/Metabolite_Brain_Gut_MorgunLab/network_properties/median_norm_data/updated_metaanalysis_fdr_table.csv")
df_q = pd.read_csv("/Users/laurenbell/Desktop/metabolomics_data/FECI/feci_analysis/feci_fdr_table.csv")


#thresholds to iterate over
ind_p = [1.0, 0.3, 0.25, 0.2, 0.15, 0.1, 0.05, 0.01]
fdr_p =[1.0, 0.3, 0.25, 0.2, 0.15, 0.1, 0.05, 0.01]

#table_making(df_quant, "quantile")
#table_making(df_median_with_zeros, "median_with_zeros")
table_making(df_q, "feci")

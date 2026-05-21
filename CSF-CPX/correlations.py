import pandas as pd
import numpy as np
import itertools
import re
from scipy import stats
import chime


def compute_correlations(data, model_name):
    nodes = data['ID'].values

    gene_data = [(i, node_val) for i, node_val in enumerate(nodes)
                 if str(node_val).startswith('ENSMUSG') and str(node_val).endswith('-P')]
    metabolite_data = [(i, node_val) for i, node_val in enumerate(nodes)
                       if not str(node_val).startswith('ENSMUSG')]

    mice = [col for col in data.columns if col != 'ID']
    pairs = itertools.product(metabolite_data, gene_data)

    results = []
    for (p_idx, p_id), (c_idx, c_id) in pairs:
        p_values = data.iloc[p_idx][mice].astype(float).to_numpy()
        c_values = data.iloc[c_idx][mice].astype(float).to_numpy()

        mask = ~np.isnan(p_values) & ~np.isnan(c_values)
        p_masked = p_values[mask]
        c_masked = c_values[mask]
        n = int(mask.sum())

        if n < 3 or np.all(p_masked == p_masked[0]) or np.all(c_masked == c_masked[0]):
            corr, pval = np.nan, np.nan
        else:
            corr, pval = stats.spearmanr(p_masked, c_masked, nan_policy='omit')
            if abs(corr) == 1:
                pval = 0

        results.append({
            'cpx_gene': c_id,
            'pls_metabolite': p_id,
            f'{model_name}_r': corr,
            f'{model_name}_p': pval,
            f'{model_name}_n': n,
        })

    return pd.DataFrame(results)


def extract_sara_id(col_name, valid_sara_ids):
    tokens = re.split(r'[\s_+]+', str(col_name))
    for token in tokens:
        token = token.strip(' ,;')
        if not token:
            continue
        if token in valid_sara_ids:
            return token
        if token.startswith('B') and token[1:] in valid_sara_ids:
            return token[1:]
    return None


def load_model_data(model_name, mouse_map):
    csf_path = f'data/CSF_{model_name}.csv'
    cpx_path = f'data/CPX_{model_name}.csv'

    csf = pd.read_csv(csf_path, index_col=0)
    csf.index.name = 'ID'
    csf.reset_index(inplace=True)

    cpx = pd.read_csv(cpx_path, index_col=0)
    cpx.index.name = 'ID'
    cpx.reset_index(inplace=True)

    cpx.columns = cpx.columns.astype(str)
    csf.columns = [str(col) for col in csf.columns]

    sara_to_matt = dict(zip(mouse_map['sara_id'].astype(str), mouse_map['matt_id'].astype(str)))
    valid_sara = set(sara_to_matt.keys())

    rename_map = {}
    for col in csf.columns:
        if col == 'ID':
            continue
        sara_id = extract_sara_id(col, valid_sara)
        if sara_id:
            rename_map[col] = sara_to_matt[sara_id]

    csf = csf.rename(columns=rename_map)

    common_cols = [col for col in csf.columns if col != 'ID' and col in cpx.columns]

    merged = pd.concat([
        cpx[['ID'] + common_cols],
        csf[['ID'] + common_cols]
    ], axis=0, ignore_index=True)

    merged.to_csv(f'merged_{model_name.lower()}.csv', index=False)
    print(f'merged {model_name}: {len(common_cols)} samples')
    return merged


def main():
    mouse_map = pd.read_csv('data/mouse_map.csv', dtype=str)
    mouse_map.columns = ['matt_id', 'sara_id']

    merged_dss = load_model_data('DSS', mouse_map)
    dss_correlations = compute_correlations(merged_dss, 'DSS')
    dss_correlations.to_csv('DSS_correlations.csv', index=False)
    print('DSS correlations done')
    chime.success()

    merged_vecpac = load_model_data('VECPAC', mouse_map)
    vecpac_correlations = compute_correlations(merged_vecpac, 'VECPAC')
    vecpac_correlations.to_csv('VECPAC_correlations.csv', index=False)
    print('VECPAC correlations done')
    chime.success()

    merged_lps = load_model_data('LPS', mouse_map)
    lps_correlations = compute_correlations(merged_lps, 'LPS')
    lps_correlations.to_csv('LPS_correlations.csv', index=False)
    print('LPS correlations done')
    chime.success()

    dss_renamed = merged_dss.rename(columns={col: f'{col}_DSS' if col != 'ID' else col for col in merged_dss.columns})
    vecpac_renamed = merged_vecpac.rename(columns={col: f'{col}_VECPAC' if col != 'ID' else col for col in merged_vecpac.columns})
    lps_renamed = merged_lps.rename(columns={col: f'{col}_LPS' if col != 'ID' else col for col in merged_lps.columns})

    vecpac_renamed.drop('ID', axis=1, inplace=True)
    lps_renamed.drop('ID', axis=1, inplace=True)

    all_data = pd.concat([dss_renamed, vecpac_renamed, lps_renamed], axis=1)
    all_data.to_csv('merged_all_data.csv', index=False)

    pooled_correlations = compute_correlations(all_data, 'pooled')
    pooled_correlations.to_csv('pooled_correlations.csv', index=False)
    print('pooled correlations done')
    chime.success()


if __name__ == '__main__':
    main()
python calc_bibc.py --network intables/FEC-STM_via_pls_edges.csv --type_map intables/FEC-STM_via_pls_typemap.csv --type1 "Feces" --type2 "Striatum" --output outtables/FEC-STM_via_pls_bibc.csv &

python calc_bibc.py --network intables/PLS-STM_edges.csv --type_map intables/PLS-STM_typemap.csv --type1 "Plasma" --type2 "Striatum" --output outtables/PLS-STM_bibc.csv &

python calc_bibc.py --network intables/FEC-STM_via_cpx_edges.csv --type_map intables/FEC-STM_via_cpx_typemap.csv --type1 "Feces" --type2 "Striatum" --output outtables/FEC-STM_via_cpx_bibc.csv &

python calc_bibc.py --network intables/FEC-STM_via_pls_cpx_edges.csv --type_map intables/FEC-STM_via_pls_cpx_typemap.csv --type1 "Feces" --type2 "Striatum" --output outtables/FEC-STM_via_both_bibc.csv &
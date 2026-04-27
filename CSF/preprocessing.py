"""
STEP 1: 
     - Raw data in Box, PROJECTS.Brain-gut_SaraCarloni.Metabolomics
     - Drop samples according to QC
     - Filter for CTRLs
     - Log transformation
     - Quantile normalization
     - Drop 0s for NaNs
"""

import pandas as pd

qc = pd.read_csv("data/CSF_QC.csv")
df = pd.read_csv("data/CSF_raw.csv", header=None)

#rename the sample names
df.iloc[0] = df.iloc[0].astype(str).str.strip().str.replace(r'  ', '_', regex=True).replace(r' ', '_', regex=True)
qc["SampleID"] = qc["SampleID"].str.strip().str.replace(r' ', '_', regex=True)
df.columns = df.iloc[0]

#if action == "excluded", drop sample
df = df.drop(columns=['CSF_101_102_104_21122021','CSF_105_21122021','CSF_106_20122021',
                      'CSF_75_103_20122021','DSS_T1.6+901_090224','DSS_BT1.8_090224',
                      'DSS_T3.7_090224','952+953_CSF_2,5H_080224','944+946_CSF_2,5H_080224'])



#if pooled == true, drop sample
#to_drop = qc[qc["Pooled"] == "TRUE"]["SampleID"] #doesnt work, figure out why 
pooled = ["DSS_71+77+T0.5_090224", "DSS_72+74+T0.8.2_090224", "DSS_T0.5+B897_090224",
"DSS_79+111+80+902_090224", "DSS_T1.5+T1.7+902_090224", "CSF_115_116_121_122_21122021",
"CSF_88_118_21122021", "DSS_122+90+87_120224", "DSS_906+T3.8_120224",
"CSF_123_124_126_20122021", "DSS_911+T5.6_120224", "DSS_93+126+T5.6_120224",
"DSS_T5.5+T5.8_120224","DSS_129+94+128_120224","DSS_95+127_120224",
"945+949_CTR_2,5H_CSF_070224", "950+951+947_CSF_2,5H_080224", "954+961_CSF_24H_080224",
"955+961_CSF_24H_080224", "960+962+963_CSF_24H_080224", "958+963_CSF_24H_080224",
"938+942_CTR_LPS_CSF_070224", "939+940.2_CTR_LPS_CSF_070224", "940+934_CTR_LPS_CSF_070224",
"3296+3296.1_VECPAC_200224", "3306+3296+3296.1VECPAC_200224_2", "3328+3356_VECPAC_190224",
"3330+3347+3689_VECPAC_190224", "3352+3370+3689_VECPAC_200224", "3355+3364+3665_VECPAC_200224",
"VECPAC_3678+3675_120224", "3369+3673+3371+3677_VECPAC_190224", "3354+3353VECPAC_190224",
"3674+3677_VECPAC_190224", "3677+3688_VECPAC_190224", "3680+3629.1_REP_VECPAC_200224",
"CSF_81_83_108_113_21122021", "DSS_114+107_090224", "CSF_118_120_21122021",
"DSS_85+86+T3.8+78+906_120224"] #you dumbass you could've done this programatically 

qced = df.drop(columns=pooled)
qced = qced[1:]
qced.to_csv("QC_filtering.csv", index=False)

#filter for controls 
#how could you do this with indicies
controls = qced.iloc[1].str.contains('CTR|T0', na=False)
columns_to_keep = controls | controls.index == 0
ctrls = qced.loc[:, columns_to_keep]
ctrls.to_csv("CTRLs_filtering.csv", index=False)


#========================================================================================

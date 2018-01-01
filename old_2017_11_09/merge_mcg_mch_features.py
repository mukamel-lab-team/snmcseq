#!/usr/bin/env python3

import pandas as pd

out_fname = './preprocessed/human_hv1_hv2_CGCH_genebody_nmcc.mat'
df_mcg = pd.read_table('./preprocessed/human_hv1_hv2_CG_genebody_nmcc.mat')
df_mch = pd.read_table('./preprocessed/human_hv1_hv2_CH_genebody_nmcc.mat')

df_merged = df_mcg.append(df_mch, ignore_index=True)
print(df_merged.shape)
print(df_mcg.shape)
print(df_mch.shape)

df_merged.to_csv(out_fname, sep='\t', header=True, index=False, na_rep='NA')

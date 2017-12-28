#!/usr/bin/env python3

import pandas as pd

input_1 = './data/binc/binc_mCH_human_combined_100000_summary_normalized.tsv'
input_2 = './data/binc/binc_mCG_human_combined_100000_summary_normalized.tsv'
output = './data/binc/binc_mCHmCG_human_combined_100000_summary_normalized.tsv'

df1 = pd.read_table(input_1)
df2 = pd.read_table(input_2)

print(df1.shape)
print(df2.shape)

df_concat = pd.concat([df1, df2], ignore_index=True)
print(df_concat.shape)

df_concat.to_csv(output, sep='\t', na_rep='NA', header=True, index=False)
print('Saved to %s' % output)
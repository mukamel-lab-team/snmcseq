#!/usr/bin/env python3

import pandas as pd
import numpy as np

df_new = pd.DataFrame()

df = pd.read_table('test_tsne_out.tsv', header=0)

for i, row in df.iterrows():
	df_new[row['cells']+'_mcc'] = np.array([row['tsne_x'], row['tsne_y']])

df_new.to_csv('test_backspin.tsv', sep='\t', header=True, index=False)
print(df_new)
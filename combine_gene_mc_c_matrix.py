#!/usr/bin/env python
"""
generate a combined matrix from matrix of each biosample
"""
import pandas as pd
import os

##### specify input #####

path = './data/gene_level'
fnames = ['genebody_mCH_MB_v1_summary.tsv', 
		'genebody_mCH_MB_EA_summary.tsv',
		'genebody_mCH_MB_EB_summary.tsv',
		]
output_fname = os.path.join(path, 'genebody_mCH_human_combined_summary.tsv') 

##### code #####

df_combined = pd.DataFrame()
cols = []
for i, fname in enumerate(fnames):
	df = pd.read_table(os.path.join(path, fname), index_col='id')
	print(df.shape)

	if i == 0:
		genes = df.index.tolist()
	else:
		assert genes == df.index.tolist()

	for col in df.columns:
		df_combined[col] = df[col]	
		cols.append(col)

# index (gene id) is automatically assigned by pandas
df_combined = df_combined[cols]
print(df_combined.shape)
df_combined.to_csv(output_fname, sep='\t', na_rep='NA', header=True, index=True)


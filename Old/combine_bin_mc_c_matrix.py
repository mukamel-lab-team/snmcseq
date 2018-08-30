#!/usr/bin/env python
"""
generate a combined matrix from matrix of each biosample
"""
import pandas as pd
import os

##### specify input #####

path = './data/binc'
fnames = ['binc_mCG_MB_v1_100000_summary.tsv', 
		'binc_mCG_MB_EA_100000_summary.tsv',
		'binc_mCG_MB_EB_100000_summary.tsv',
		]
output_fname = os.path.join(path, 'binc_mCG_human_combined_100000_summary.tsv') 

##### code #####

df_combined = pd.DataFrame()
cols = []
for i, fname in enumerate(fnames):
	df = pd.read_table(os.path.join(path, fname), 
		index_col=['chr', 'bin'], 
		dtype={'chr': object, 'bin': int})
	print(df.shape)

	if i == 0:
		chr_bin = df.index.tolist()
	else:
		assert chr_bin == df.index.tolist()

	for col in df.columns:
		df_combined[col] = df[col]	
		cols.append(col)

# index is automatically assigned by pandas
df_combined = df_combined[cols]
print(df_combined.shape)
df_combined.to_csv(output_fname, sep='\t', na_rep='NA', header=True, index=True)
print('Saved to %s' % output_fname)


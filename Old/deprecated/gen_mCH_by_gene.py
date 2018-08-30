#!/usr/bin/env python3

import pandas as pd
import numpy as np
import glob
import os 

DIR = '/cndd/Public_Datasets/single_cell_methylome/gene_level/human/genebody_MB_EB'
metadata = '/cndd/fangming/side_project/snmcseq_dev/metadata/MB_EB_metadata_cells.tsv'
# test_samples = ['Pool_2256_AD002_indexed_R1', 
# 				'Pool_2256_AD006_indexed_R1',
# 				'Pool_2256_AD008_indexed_R1', 
# 				'Pool_2256_AD010_indexed_R1'] 

df_meta = pd.read_table(metadata, header=0)
df_meta = df_meta[['Sample', 'mCH/CH']]
# df_meta = df_meta[df_meta['Sample'].isin(test_samples)]
df_meta.sort_values('Sample', axis=0, inplace=True)

# read all in 
df_combined = pd.DataFrame()
fnames = glob.glob(os.path.join(DIR, '*.txt'))
for fname in fnames:
	print('.', end='', flush=True)
	sample = os.path.basename(fname)
	if sample.endswith('_mch_genebody.txt'):
		sample = sample[:-len('_mch_genebody.txt')]
	df = pd.read_table(fname, header=0)		
	df['cells'] = np.array([sample]*df.shape[0])
	df_combined = pd.concat([df_combined, df])

# loop through all genes
# generate files
for i, group in enumerate(df_combined.groupby('id')):
	# format: id/cell/mCH/normalized_mCH (automatic generate 'nan' for 0/0 case)
	print('*', end='', flush=True)
	df_new = pd.DataFrame()
	label, df = group
	df.sort_values('cells', axis=0, inplace=True)
	assert df['cells'].tolist() == df_meta['Sample'].tolist()

	df_new['id'] = df['id']
	df_new['cells'] = df['cells']
	df_new['mcc'] = df['mc'].values/df['c'].values
	df_new['norm_mcc'] = df_new['mcc'].values/df_meta['mCH/CH'].values
	df_new = df_new[['id', 'cells', 'mcc', 'norm_mcc']]
	df_new.to_csv('./data/%s_mCH.txt' % label, sep='\t', header=False, index=False)


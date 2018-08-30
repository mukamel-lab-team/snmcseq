#!/usr/bin/env python3

import pandas as pd
import numpy as np
import glob
import os 

DIR = '/cndd/Public_Datasets/single_cell_methylome/gene_level/human/genebody_v1'
metadata = '/cndd/fangming/side_project/snmcseq_dev/metadata/human_v1_metadata_cells.tsv'
OUT_DIR = '/cndd/Public_Datasets/single_cell_methylome/gene_level/human/genebody_v1_by_gene'
 
df_meta = pd.read_table(metadata, header=0)
df_meta = df_meta[['Sample', 'mCH/CH']]
### for test purpose
# test_samples = ['Pool_2256_AD002_indexed_R1', 
# 				'Pool_2256_AD006_indexed_R1',
# 				'Pool_2256_AD008_indexed_R1', 
# 				'Pool_2256_AD010_indexed_R1'] 
# df_meta = df_meta[df_meta['Sample'].isin(test_samples)]

# sort in the order of sample name
df_meta.sort_values('Sample', axis=0, inplace=True)
df_meta = df_meta.reset_index(drop=True)

# read all in 
dfs = []
fnames = glob.glob(os.path.join(DIR, '*.txt'))
# read in the order of sample name
fnames = sorted(fnames)
for i, fname in enumerate(fnames):
	print('.', end='', flush=True)
	sample = os.path.basename(fname)
	if sample.endswith('_mch_genebody.txt'):
		sample = sample[:-len('_mch_genebody.txt')]
	df = pd.read_table(fname, header=0)		
	df.sort_values('id', inplace=True) ### very important!!!
	df = df.reset_index(drop=True) ### equally important!!!
	# do check here
	# if i == 0:
	# 	gene_ids = df['id'].tolist()		
	# else:
	# 	assert df['id'].tolist() == gene_ids
	df['cells'] = np.array([sample]*df.shape[0])
	df['mcc'] = df['mc'].values/df['c'].values
	df = df[['id', 'cells', 'mcc']]
	dfs.append(df)

# nrow: numher of genes
# ncol: number of attributes 
# each df is a sample/cell
nrow, ncol = dfs[0].shape
# check dimension
for df in dfs:
	assert df.shape == (nrow, ncol)	
# check sample order
sample_order = []
for df in dfs:
	sample_order.append(df.ix[0,'cells'])
print(len(sample_order))
print(len(df_meta['Sample'].tolist()))
print(sample_order[:100])
print(df_meta['Sample'].tolist()[:100])
assert sample_order == df_meta['Sample'].tolist()
# check gene order
for i, df in enumerate(dfs):
	if i == 0:
		gene_ids = df['id'].tolist()
	else:
		assert df['id'].tolist() == gene_ids


# loop through all genes
# i is gene; j is sample
#for i in range(10):
for i in range(nrow):
	print('*', end='', flush=True)
	dict_list = []
	gene_id = dfs[0].ix[i,'id']
	# check gene order   ### very important !!!
	for df in dfs:
		assert df.ix[i,'id'] == gene_id
	# loop through all samples
	for j, df in enumerate(dfs):
		mcc = df.ix[i,'mcc']
		mcc_norm = mcc/df_meta.ix[j,'mCH/CH'] 
		dict_info = {'id': gene_id, 
					'cells': df.ix[i,'cells'],
					'mcc': mcc, 
					'mcc_norm': mcc_norm}
		dict_list.append(dict_info)	
	df_gene = pd.DataFrame(dict_list)
	df_gene = df_gene[['id', 'cells', 'mcc', 'mcc_norm']]
	df_gene.to_csv(os.path.join(OUT_DIR, '%s_mCH.txt' % gene_id),
		sep='\t', na_rep='nan', header=False, index=False)

print('Done!')

#!/usr/bin/env python3

import pandas as pd
import numpy as np
import glob
import os 

metadata = '/cndd/fangming/snmcseq_dev/data/metadata/metadata_human_combined_updated.tsv'
input_matrix = '/cndd/fangming/snmcseq_dev/data/gene_level/genebody_mCG_human_combined_summary.tsv'
output_dir = '/cndd/fangming/snmcseq_dev/data/gene_level/genebody_mCG_by_gene_human_combined'
context = 'CG'




# global mcc from metadata 
print('Processing meta data...')
if context=='CG':
	sr_global_mcc = pd.read_table(metadata, header=0, index_col='Sample')['mCG/CG'] 
elif context=='CH':
	sr_global_mcc = pd.read_table(metadata, header=0, index_col='Sample')['mCH/CH'] 

sr_global_mcc.index = [index+'_mcc' for index in sr_global_mcc.index]

# mc_c mtrx
print('Loading data...')
# for test purpose
# df_mc_c_mtrx = pd.read_table(input_matrix, header=0, usecols=[i for i in range(10)], index_col='id')
# sr_global_mcc = sr_global_mcc.loc[[col[:-len('_c')]+'_mcc' for col in df_mc_c_mtrx.filter(regex='_c$').columns]]
# for real run
df_mc_c_mtrx = pd.read_table(input_matrix, header=0, index_col='id')


print('Computing mcc...')
# mcc mtrx
df_mc = df_mc_c_mtrx.filter(regex='_mc$')
df_mc.columns = [col[:-len('_mc')]+'_mcc' for col in df_mc.columns]
df_c = df_mc_c_mtrx.filter(regex='_c$')
df_c.columns = [col[:-len('_c')]+'_mcc' for col in df_c.columns]
df_mcc = df_mc/df_c 

print('Computing normalized mcc...')
# nmcc mtrx
df_norm_mcc = df_mcc.divide((sr_global_mcc+0.0001), axis=1)
df_norm_mcc.columns = [col[:-len('_mcc')]+'_nmcc' for col in df_norm_mcc.columns]

# check if both mtrx has the same index and columns
assert ([col[:-len('_mcc')] for col in df_mcc.columns] == 
	[col[:-len('_nmcc')] for col in df_norm_mcc.columns])
assert df_mcc.index.tolist() == df_norm_mcc.index.tolist()

print('Outputing...')
# output 
for gene_id, row_mcc in df_mcc.iterrows():
	df_gene = pd.DataFrame(columns=['gene_id', 'cell', 'mcc', 'norm_mcc']) 
	df_gene['cell'] = [col[:-len('_mcc')] for col in df_mcc.columns]
	df_gene['mcc'] = row_mcc.values
	df_gene['norm_mcc'] = df_norm_mcc.loc[gene_id, :].values 
	df_gene['gene_id'] = gene_id

	if not os.path.exists(output_dir):
	    os.makedirs(output_dir)

	output_fname = os.path.join(output_dir, '%s_m%s.txt' % (gene_id, context)) 
	df_gene.to_csv(output_fname, 
			sep='\t', na_rep='NA', index=False, header=False)
	print(df_gene.shape)
	print('Saved to %s.' % output_fname)

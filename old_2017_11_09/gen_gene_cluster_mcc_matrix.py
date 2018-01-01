#!/usr/bin/env python3

"""
generate gene_cluster_mcc_matrix from gene_cell_mc_c_matrix


input: gene_cell_mc_c_matrix, cell_metadata

metadata is used to get cluster information

INDEPENDENT OF mCG or mCH

output: gene_cluster_mcc_matrix
"""

import pandas as pd
import os
import glob
import my_utils

input_fname = './preprocessed/genebody_mCG_combined-hs1-hs2.mat'
meta_fname = './metadata/human_hv1_hv2_metadata_cells.tsv'
base_call_cutoff = 50
output_fname = './preprocessed/glia/human_hv1_hv2_gene_cluster_vCG_mcg_matrix.tsv'

input_cluster = './preprocessed/human_hv1_hv2_CG_genebody_nmcc_backspin.tsv'
NEW_BACKSPIN = True

# logger
logger = my_utils.create_logger() 

# get cluster_cell info from meta
logger.info('Loading data...')
df_meta = pd.read_table(meta_fname, header=0, index_col=None)
df_meta = df_meta[['Sample', 'cluster_ID']]
if NEW_BACKSPIN:
	df_meta = df_meta.drop('cluster_ID', axis=1)
	df_cluster = pd.read_table(input_cluster, header=0)
	df_meta = pd.merge(df_meta, df_cluster, left_on='Sample', right_on='Sample', how='inner')

dict_cluster = {}
for cluster_ID, df_meta_sub in df_meta.groupby('cluster_ID'):
	dict_cluster[cluster_ID] = df_meta_sub.Sample.tolist()
cluster_ids = list(dict_cluster.keys())
cluster_ids = sorted([int(item[len('cluster_'):]) for item in cluster_ids])
cluster_ids = ['cluster_'+str(item) for item in cluster_ids]

# get info from df
# for test purpose
# df = pd.read_table(input_fname, header=0, index_col=0, nrows=100)
# for real run
df = pd.read_table(input_fname, header=0, index_col=0)

# get samples 
samples = [col[:-len('_c')] for col in df.columns if col.endswith('_c')]

# filter out samples not in metadata 
samples = [sample for sample in samples if sample in df_meta.Sample.tolist()]
df_cols = []
for sample in samples:
	df_cols.append(sample+'_mc')
	df_cols.append(sample+'_c')
df = df[df_cols]	

### filter out genes
logger.info('filter out genes and calculate...')
# generate 2 information matrices: 
# coverage matrix is used to generate information matrix
# cluster_ID is NOT in order!
# coverage matrix
df_coverage = pd.DataFrame(index=df.index)
for cluster_ID, cluster_samples in dict_cluster.items():

	cluster_sample_mc_cols = [sample+'_mc' for sample in cluster_samples] 
	cluster_sample_c_cols = [sample+'_c' for sample in cluster_samples] 
	# coverage matrix
	cluster_size = len(cluster_samples)
	df_coverage[cluster_ID+'_nc'] = (
		(df.filter(regex='_c$')[cluster_sample_c_cols]>base_call_cutoff).sum(axis=1)
		)
	df_coverage[cluster_ID+'_rc'] = df_coverage[cluster_ID+'_nc']/cluster_size 

condition1 = df_coverage.filter(regex='_nc$').min(axis=1) > 10  # all cell > 10
condition2 = df_coverage.filter(regex='_rc$').min(axis=1) > 0.2 # all cell > 20%
condition3 = df_coverage.filter(regex='_rc$').max(axis=1) > 0.5 # at least 1 cluster > 50%

logger.info('gene*cells before pre-filtering: ' + str(df.shape))
df = df[(condition1 & condition2 & condition3)]
logger.info('gene*cells after pre-filtering: ' + str(df.shape))

### doing mcc calculation after filtering out genes
# information matrix
df_mcc = pd.DataFrame(index=df.index)
# mean mccs averaged across all cells 
mean_mccs = df.filter(regex='_mc$').sum(axis=1)/df.filter(regex='_c$').sum(axis=1)
df_mcc['mean_mcc'] = mean_mccs
for cluster_ID, cluster_samples in dict_cluster.items():

	cluster_sample_mc_cols = [sample+'_mc' for sample in cluster_samples] 
	cluster_sample_c_cols = [sample+'_c' for sample in cluster_samples] 
	# information matrix
	df_mcc[cluster_ID+'_mcc'] = (
		df.filter(regex='_mc$')[cluster_sample_mc_cols].sum(axis=1)/
		df.filter(regex='_c$')[cluster_sample_c_cols].sum(axis=1)
		)

# # normalize across all genes
# df_mcc = df_mcc.divide(mean_mccs, axis=0)
# reorder columns
df_mcc = df_mcc[['mean_mcc']+[cluster_ID+'_mcc' for cluster_ID in cluster_ids]]
logger.info('output gene*cluster matrix: ' + str(df_mcc.shape))
print(df_mcc.head())

# output to file
logger.info('Saving to file...')
df_mcc.to_csv(output_fname, 
	sep='\t', na_rep='NA', header=True, index=True)
logger.info('Saved to %s' % output_fname)


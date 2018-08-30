#!/usr/bin/env python3
"""
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import seaborn as sns
from scipy.cluster import hierarchy

def gene_name_to_id(name, df_gene_id):
	try:
		gene_id = df_gene_id[df_gene_id.geneName==name].geneID.values[0]
		return gene_id
	except:
		print('Gene %s not found!' % name)
		# raise ValueError('Gene %s not found!' % name)

def gene_id_to_name(gene_id, df_gene_id):
	try:
		gene_name = df_gene_id[df_gene_id.geneID==gene_id].geneName.values[0]
		return gene_name
	except:
		print('Gene %s not found!' % gene_id)
		# raise ValueError('Gene %s not found!' % name)

def np_mean_with_percentiles(array, low_p=5, high_p=95):
	"""
	"""
	lo = np.nanpercentile(array, low_p)
	hi = np.nanpercentile(array, high_p)
	array = array[array > lo]
	array = array[array < hi]

	return np.nanmean(array)
	

def get_hc_clusters(array):
	"""
	2-column step by step clustering info

	array should be int ndarray
	"""

	nrow, ncol = array.shape
	assert ncol == 2	

	# a list where each entry is a list of cluster numbers (start from 1)
	cluster_info = [[i] for i in range(1, nrow+1+1)]

	for row in array:
		# array matrix gives info about cluster indices
		# a_ind is where a_ind+1 cluster is stored
		a_ind, b_ind = row[0], row[1]

		a_cluster = cluster_info[a_ind]
		b_cluster = cluster_info[b_ind]
		new_cluster = sorted(a_cluster + b_cluster)
		cluster_info.append(new_cluster)

	return cluster_info


# look at a few marker genes
marker_genes = ['Grm3', 'Slc6a11', 'Fgfr3', 'Slc4a4', 'Mlc1',  # Astrocytes
				# 'Reln', 'Snhg11', 'Islr2', 'Npy', 'Tmem130',  # Neurons
				'Pdgfra', 'Mmp15', 'Cdo1', 'Rlbp1', 'Lhfpl3',  # OPCs
				'Fyn', 'Nfasc', 'Enpp6', 'Tmem163', 'Lims2',  # NFO
				'Gjb1', 'Ndrg1', 'Ppp1r14a', 'Mbp', 'Mal',  # MO
				'Gpr84', 'Tnf', 'Ncf1', 'Gdf15', 'Osm',  # Microglia 
				# 'Cldn5', 'Bsg', 'Slc16a1', 'Vwa1', 'Cd34', # Endothelial
				# 'Igf2', 'Rdh10', 'Col1a1', 'Dcn', 'Col1a2'  # Pericytes
				]
marker_genes = [gene.upper() for gene in marker_genes]

gene_id_fname = './metadata/gene_id_to_names.tsv'
df_gene_id = pd.read_table(gene_id_fname, header=0)

marker_gene_ids = [gene_name_to_id(gene, df_gene_id)
				for gene in marker_genes]


# glia_clusters from metadata
input_meta = './metadata/human_hv1_hv2_metadata_cells.tsv' 
df_meta = pd.read_table(input_meta, header=0)
df_meta = df_meta[['Sample', 'cluster_ID', 'cluster_label']]
df_meta = df_meta.loc[df_meta.cluster_label=='glia', :]
# glia cluster_ID list
cluster_IDs = np.unique(df_meta['cluster_ID'].values)

# get gene*cluster normalized mcc (across all cells for a gene)
input_mcc = './preprocessed/human_hv1_hv2_gene_cluster_mch_matrix.tsv'
df_mcc = pd.read_table(input_mcc, header=0, index_col='id')
df_mcc = df_mcc[[cluster_ID+'_mcc' for cluster_ID in cluster_IDs]]
df_mcc = df_mcc.loc[df_mcc.index.isin(marker_gene_ids), :]

df_mcc['geneName'] = [gene_id_to_name(gene_id, df_gene_id)
				for gene_id in df_mcc.index]

df_mcc.set_index('geneName', inplace=True)
print(df_mcc.shape)


# row_linkage = hierarchy.linkage(df_mcc.values, method='average')
# col_linkage = hierarchy.linkage(df_mcc.values.T, method='average')
# # dendro = hierarchy.dendrogram(col_linkage)

g = sns.clustermap(df_mcc, 
				col_cluster=True, row_cluster=True, 
				# row_linkage = row_linkage, col_linkage = col_linkage,
				figsize=(8,6), cmap='viridis')
plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)
plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90)
g.savefig('./preprocessed/cluster_glias.pdf')

plt.show()

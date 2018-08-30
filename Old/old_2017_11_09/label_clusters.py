#!/usr/bin/env python3
"""
Give each cluster a label: 
[excitatory, inhibitory, glia, unknown]
based on marker gene normalized mCH level
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import seaborn as sns
from scipy.spatial import distance
from scipy.cluster import hierarchy

def gene_name_to_id(name, df_gene_id):
	try:
		gene_id = df_gene_id[df_gene_id.geneName==name].geneID.values[0]
		return gene_id
	except:
		raise ValueError('No such gene found!')

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
marker_genes_exci = ['SATB2', 'TYRO3', 'ARPP21', 'SLC17A7', 'TBR1', 'CAMK2A', 'ITPKA']
marker_genes_inhi = ['GAD1', 'ERBB4', 'SLC6A1']
marker_genes_not_glia = ['MEF2C']
marker_genes = (marker_genes_exci + 
				marker_genes_inhi +
				marker_genes_not_glia)

gene_id_fname = './metadata/gene_id_to_names.tsv'
df_gene_id = pd.read_table(gene_id_fname, header=0)

marker_gene_ids_exci = [gene_name_to_id(gene, df_gene_id)
				for gene in marker_genes_exci]
marker_gene_ids_inhi = [gene_name_to_id(gene, df_gene_id)
				for gene in marker_genes_inhi]
marker_gene_ids_not_glia = [gene_name_to_id(gene, df_gene_id)
				for gene in marker_genes_not_glia]
marker_gene_ids = [gene_name_to_id(gene, df_gene_id)
				for gene in marker_genes]


# # heatmap data
# heat_matrix = df[df.name.isin(marker_genes)][cluster_cols_ordered].values

### get gene*cluster info
# gene mCH files
gene_dir = ('/cndd/Public_Datasets/single_cell_methylome'
			+'/gene_level/human/genebody_combined_by_gene')
# cluster info
cluster_fname = './preprocessed/human_hv1_hv2_CH_100000_clusters.tsv'
# format: sample, cluster_ID
df_cluster = pd.read_table(cluster_fname, header=0)

# loop over all selected genes
info_dict_list = []
for (marker_gene_id, marker_gene) in zip(marker_gene_ids, marker_genes):
	# format: no header; gene_id/sample/mcc/norm_mcc
	df_gene = pd.read_table(os.path.join(gene_dir, marker_gene_id+'_mCH.txt'), header=None)
	df_gene.columns = ['gene_id', 'sample', 'mcc', 'norm_mcc']

	df_gene_cluster = pd.merge(df_gene, df_cluster, left_on='sample', right_on='sample')
	df_gene_cluster = df_gene_cluster[['cluster_ID', 'sample', 'norm_mcc']]
	# loop over all clusters
	for cluster_ID, df_gene_cluster_sub in df_gene_cluster.groupby('cluster_ID'):
		# marker_gene, cluster, norm_mccs
		norm_mccs = df_gene_cluster_sub.norm_mcc.values
		mean_mcc = np_mean_with_percentiles(norm_mccs)
		info_dict_list.append({
							'gene_ID': marker_gene, 
							'cluster_ID': cluster_ID[len('cluster_'):], 
							'mean_mcc': mean_mcc})

df_info = pd.DataFrame(info_dict_list)
df_info = pd.pivot_table(df_info, values='mean_mcc', index=['gene_ID'], columns=['cluster_ID'])

# re-order the index of df columns (very important)!
print(df_info.columns)
cluster_cols = df_info.columns.tolist()
cluster_nums = sorted([int(item) 
 			for item in cluster_cols])
cluster_cols_ordered = ['%d' % item
			for item in cluster_nums]
df_info = df_info[cluster_cols_ordered]
print(df_info.columns)


row_linkage = hierarchy.linkage(df_info.values, method='average')
col_linkage = hierarchy.linkage(df_info.values.T, method='average')
# dendro = hierarchy.dendrogram(col_linkage)

# get clustering relationship 
n_genes, n_clusters = df_info.shape
hc_info = col_linkage.astype('int')[:,:2]
hc_results_info = get_hc_clusters(hc_info)
# print([len(item) for item in hc_results_info])

# get cluster numbers (hard coded)
for hc_result in hc_results_info:
	if len(hc_result) == 13: # glia
		cluster_a = hc_result
		print('found cluster a')
	elif len(hc_result) == 28: # exci
		cluster_b = hc_result
		print('found cluster b')
	elif len(hc_result) == 41: # inhi
		cluster_c = hc_result
		print('found cluster c')

# have to manually check if we get the right cluster labels
print(cluster_a) # glia
print(cluster_b) # exci
print(cluster_c) # inhi

# output to file
output_fname_hc = './preprocessed/human_hv1_hv2_CH_100000_cluster_labels.tsv'
dict_list = []
for cluster_num in cluster_a:
	dict_list.append({'cluster_num': cluster_num, 'cluster_label': 'glia'})
for cluster_num in cluster_b:
	dict_list.append({'cluster_num': cluster_num, 'cluster_label': 'exci'})
for cluster_num in cluster_c:
	dict_list.append({'cluster_num': cluster_num, 'cluster_label': 'inhi'})
df = pd.DataFrame(dict_list)
df.sort_values(['cluster_num'], inplace=True)
df.to_csv(output_fname_hc, sep='\t', header=True, index=False)


# g = sns.clustermap(df_info, col_cluster=True, row_cluster=True, 
# 				row_linkage = row_linkage, col_linkage = col_linkage,
# 				figsize=(8,6), cmap='viridis')
# plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)
# plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90)
# g.savefig('./preprocessed/cluster_exci_inhi_glia.pdf')

plt.show()



















































###### deprecated #######
# cluster_cols = df_info.columns.tolist()
# cluster_nums = sorted([int(item) 
#  			for item in cluster_cols])
# cluster_cols_ordered = ['%d' % item
# 			for item in cluster_nums]
# df_info = df_info[cluster_cols_ordered]


# df_info.ix[:, (df_info.loc['MEF2C', :] > 0.5)] = np.nan

# df_exci = df_info[df_info.index.isin(marker_genes_exci)]
# df_inhi = df_info[df_info.index.isin(marker_genes_inhi)]
# df_not_glia = df_info[df_info.index.isin(marker_genes_not_glia)]
# marker_genes_ordered = df_exci.index.tolist() + df_inhi.index.tolist() + df_not_glia.index.tolist() 

# heat_matrix = np.concatenate((df_exci.values, df_inhi.values, df_not_glia.values), axis=0) 

# # zscore
# means = np.mean(heat_matrix, axis=1)
# means = np.repeat(np.reshape(means, (len(means), 1)), len(cluster_cols), axis=1)
# stds = np.std(heat_matrix, axis= 1) 
# stds = np.repeat(np.reshape(stds, (len(stds), 1)), len(cluster_cols), axis=1)
# heat_matrix_zscore = (heat_matrix - means)/stds

# # heatmap
# fig, ax = plt.subplots(figsize=(10,6))
# im = ax.imshow(heat_matrix)
# fig.colorbar(im, ax=ax)
# ax.set_aspect('auto')

# ax.set_xticks(np.arange(0, len(cluster_cols_ordered)-0.5))
# xticklabels = [i if i%5==0 else '' for i in range(1, len(cluster_cols_ordered)+1)]
# ax.set_title('Normalized mCH level of marker genes (excitatory, inhibitory, and glia)')
# ax.set_xticklabels(xticklabels, rotation=90)
# ax.set_yticks(range(len(marker_genes_ordered)))
# ax.set_yticklabels(marker_genes_ordered)


# fig, ax = plt.subplots(figsize=(10,6))
# im = ax.imshow(heat_matrix_zscore)
# fig.colorbar(im, ax=ax)
# ax.set_aspect('auto')

# ax.set_xticks(np.arange(0, len(cluster_cols_ordered)-0.5))
# xticklabels = [i if i%5==0 else '' for i in range(1, len(cluster_cols_ordered)+1)]
# ax.set_title('Zscore of normalized mCH level of marker genes (excitatory, inhibitory, and glia)')
# ax.set_xticklabels(xticklabels, rotation=90)
# ax.set_yticks(range(len(marker_genes_ordered)))
# ax.set_yticklabels(marker_genes_ordered)





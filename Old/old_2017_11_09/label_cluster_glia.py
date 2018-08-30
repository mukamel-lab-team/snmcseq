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
		print('Gene %s not found!' % name)
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
		# only look at glia cluster_IDs
		if cluster_ID in cluster_IDs:
			# marker_gene, cluster, norm_mccs
			norm_mccs = df_gene_cluster_sub.norm_mcc.values
			mean_mcc = np_mean_with_percentiles(norm_mccs)
			if np.isnan(mean_mcc):
				print(marker_gene)

			info_dict_list.append({
								'gene_ID': marker_gene, 
								'cluster_ID': cluster_ID[len('cluster_'):], 
								'mean_mcc': mean_mcc})

df_info = pd.DataFrame(info_dict_list)
df_info = pd.pivot_table(df_info, values='mean_mcc', index=['gene_ID'], columns=['cluster_ID'])

# re-order the index of df columns (very important)!
cluster_cols = df_info.columns.tolist()
cluster_nums = sorted([int(item) 
 			for item in cluster_cols])
cluster_cols_ordered = ['%d' % item
			for item in cluster_nums]
df_info = df_info[cluster_cols_ordered]


row_linkage = hierarchy.linkage(df_info.values, method='average')
col_linkage = hierarchy.linkage(df_info.values.T, method='average')
# # dendro = hierarchy.dendrogram(col_linkage)

# ### get clustering relationship 
# # n_genes, n_clusters = df_info.shape
# # hc_info = col_linkage.astype('int')[:,:2]
# # hc_results_info = get_hc_clusters(hc_info)
# # # print([len(item) for item in hc_results_info])

# # # get cluster numbers (hard coded)
# # for hc_result in hc_results_info:
# # 	if len(hc_result) == 13: # glia
# # 		cluster_a = hc_result
# # 		print('found cluster a')
# # 	elif len(hc_result) == 28: # exci
# # 		cluster_b = hc_result
# # 		print('found cluster b')
# # 	elif len(hc_result) == 41: # inhi
# # 		cluster_c = hc_result
# # 		print('found cluster c')

# # # have to manually check if we get the right cluster labels
# # print(cluster_a) # glia
# # print(cluster_b) # exci
# # print(cluster_c) # inhi

# # # output to file
# # output_fname_hc = './preprocessed/human_hv1_hv2_CH_100000_cluster_labels.tsv'
# # dict_list = []
# # for cluster_num in cluster_a:
# # 	dict_list.append({'cluster_num': cluster_num, 'cluster_label': 'glia'})
# # for cluster_num in cluster_b:
# # 	dict_list.append({'cluster_num': cluster_num, 'cluster_label': 'exci'})
# # for cluster_num in cluster_c:
# # 	dict_list.append({'cluster_num': cluster_num, 'cluster_label': 'inhi'})
# # df = pd.DataFrame(dict_list)
# # df.sort_values(['cluster_num'], inplace=True)
# # df.to_csv(output_fname_hc, sep='\t', header=True, index=False)


g = sns.clustermap(df_info, 
				col_cluster=True, row_cluster=False, 
				row_linkage = row_linkage, col_linkage = col_linkage,
				figsize=(8,6), cmap='viridis')
plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)
plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90)
g.savefig('./preprocessed/cluster_glias.pdf')

plt.show()

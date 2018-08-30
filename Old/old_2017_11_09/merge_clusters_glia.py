#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import itertools
from scipy.cluster import hierarchy 
# from scipy.stats import zscore
# from scipy.stats import spearmanr
# from scipy.stats import entropy

# from snmcseq import plot_utils

NORMALIZE_across_cells = True 
GLIA_ONLY = True 
HC_MERGE = True 
MG_MERGE = True

### --- load data
# glia_clusters from metadata
input_meta = './metadata/human_hv1_hv2_metadata_cells.tsv' 
df_meta = pd.read_table(input_meta, header=0, index_col='Sample')
df_meta = df_meta[['cluster_ID', 'cluster_label', 'mCG/CG', 'mCH/CH']]
# glia cluster_ID list
cluster_IDs_glia = np.unique(df_meta.loc[df_meta.cluster_label=='glia', 'cluster_ID'].values)
cluster_IDs_glia = sorted([int(cluster_ID.strip('cluster_')) 
		for cluster_ID in cluster_IDs_glia])
cluster_IDs_glia = ['cluster_'+str(cluster_ID) 
		for cluster_ID in cluster_IDs_glia]

cluster_IDs = np.unique(df_meta.loc[:, 'cluster_ID'].values)
cluster_IDs = sorted([int(cluster_ID.strip('cluster_')) 
		for cluster_ID in cluster_IDs])
cluster_IDs = ['cluster_'+str(cluster_ID) 
		for cluster_ID in cluster_IDs]
# if only plot glial clusters
if GLIA_ONLY:
	cluster_IDs = cluster_IDs_glia
	df_meta = df_meta[df_meta.cluster_label=='glia']

# filtered gene-cluster mCH
input_mcc = './preprocessed/human_hv1_hv2_gene_cluster_mcg_matrix.tsv'
df_mcc = pd.read_table(input_mcc, header=0, index_col='id')
mean_mccs = df_mcc.mean_mcc
df_mcc.drop('mean_mcc', axis=1, inplace=True)
df_mcc = df_mcc[[item+'_mcc' for item in cluster_IDs]]
print(df_mcc.shape)
## normalize across all cells
if NORMALIZE_across_cells:
	df_mcc = df_mcc.divide(mean_mccs, axis=0)
# df_mcc = df_mcc.subtract(mean_mccs, axis=0)
# df_mcc = df_mcc.divide(df_mcc.std(axis=1), axis=0)
# df_mcc = df_mcc.apply(zscore, axis=1)
# ----

# gene id to name
input_gene_id_name = './metadata/gene_id_to_names.tsv'
df_gene_id_name = pd.read_table(input_gene_id_name, header=0, index_col='geneID')
df_gene_id_name = df_gene_id_name[['geneName']]

### end of data loading -------

### merge cluster (hierarchical clustering) --
if HC_MERGE:
	col_linkage = hierarchy.linkage(df_mcc.T, method='single')
	fig, ax = plt.subplots()
	dendro = hierarchy.dendrogram(col_linkage, 
								labels=[col[len('cluster_'):-len('_mcc')] for col in df_mcc.columns], 
								ax=ax)
	ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
	ax.set_title('Glial clusters; clustering with mCG levels of %d genes' % df_mcc.shape[0])
	ax.set_xlabel('cluster_ID')
	ax.set_ylabel('Single linkage, Euclidean distance')
	fig.tight_layout()
	fig.savefig('./preprocessed/hc_glia.pdf')
	print('Saved to ./preprocessed/hc_glia.pdf')
	plt.show()

### merge cluster procedure: marker gene ---
if MG_MERGE:
	CUTOFF = 20
	while True:
		num_sigs = []
		# do pair-wise marker genes 
		df_mcc_pct = df_mcc.rank(axis=0, pct=True) # axis=0 among gene
		for cluster_i, cluster_j in itertools.combinations(cluster_IDs, 2):
			# conditions 
			cond_1_i = (df_mcc_pct[cluster_i+'_mcc']<0.02)
			cond_1_j = (df_mcc_pct[cluster_j+'_mcc']>0.6)
			cond_2_i = (df_mcc_pct[cluster_i+'_mcc']>0.6)
			cond_2_j = (df_mcc_pct[cluster_j+'_mcc']<0.02)
			cond = (cond_1_i & cond_1_j) | (cond_2_i & cond_2_j)
			# marker genes
			df_spec_ij = df_mcc.loc[cond, [cluster_i+'_mcc', cluster_j+'_mcc']]	
			# sort according to diff	
			df_spec_ij['diff'] = df_spec_ij[cluster_i+'_mcc'] - df_spec_ij[cluster_j+'_mcc']
			df_spec_ij['diff'] = df_spec_ij['diff'].abs() 
			df_spec_ij = df_spec_ij.sort_values('diff', ascending=False)
			# top marker genes
			marker_genes_ij = df_spec_ij.head(min(CUTOFF, df_spec_ij.shape[0])).index 
			# print(len(marker_genes_ij))	

			# test significance stats
			# examine each marker genes, get mCH of each cell and do a t-test

			num_sigs.append({'clusters': (cluster_i, cluster_j), 
							'num_sigs': len(marker_genes_ij)})

		num_sigs = pd.DataFrame(num_sigs)
		min_num_sigs = num_sigs['num_sigs'].min()
		if min_num_sigs >= 10 or len(cluster_IDs)<=3: # exit if meeting exit standard
			break
		else: # merge 2 clusters with min(num_sigs)	
			# if duplicate, merge only 1 set of clusters
			cluster_i, cluster_j = num_sigs.loc[num_sigs['num_sigs']==min_num_sigs, 
									'clusters'].head(1).values[0]
			print(cluster_i, cluster_j)
			# merge clusters and update df_mcc and cluster_IDs
			# add 1 cluster
			df_mcc['%s-%s_mcc' % (cluster_i, cluster_j)] = df_mcc[[cluster_i+'_mcc', 
															cluster_j+'_mcc']].mean(axis=1) 
			# remove 2 clusters
			df_mcc = df_mcc.drop([cluster_i+'_mcc', cluster_j+'_mcc'], axis=1)
			# update cluster_IDs
			cluster_IDs = [col[:-len('_mcc')] for col in df_mcc.columns]
	# results
	print(cluster_IDs)
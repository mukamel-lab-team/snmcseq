#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from scipy.stats import zscore
from scipy.stats import pearsonr
from scipy.stats import spearmanr

OPC = 'Oligodendrocyte Precursor Cell'
NFO = 'Newly Formed Oligodendrocyte'
MO = 'Myelinating Oligodendrocytes'
EC = 'Endothelial Cells'
num_cutoff = 30

def rename_celltypes(cell_types):
	"""
	take a list of cell types and symplify naming convention
	"""
	new_cell_types = []
	for cell_type in cell_types:
		if cell_type == OPC:
			cell_type = 'OPC'
		elif cell_type == NFO:
			cell_type = 'NFO'
		elif cell_type == MO:
			cell_type = 'MO'
		else:
			pass
		new_cell_types.append(cell_type)
	return new_cell_types


### load data
# gene-expression fpkm
input_rna = './metadata/barreslab_rnaseq.tsv'
df_rna = pd.read_table(input_rna, header=0, index_col='Gene symbol')
# apply upper case to all entries of gene symbols
df_rna.index = df_rna.index.str.upper() 
# cell types
cell_types = df_rna.columns.tolist()
cell_types.remove('Description')
df_rna = df_rna[cell_types]
df_rna.columns = rename_celltypes(cell_types) 
# take log10
df_rna = df_rna.applymap(np.log10)
print(cell_types)
print(df_rna.head())

# glia_clusters from metadata
input_meta = './metadata/human_hv1_hv2_metadata_cells.tsv' 
df_meta = pd.read_table(input_meta, header=0)
df_meta = df_meta[['Sample', 'cluster_ID', 'cluster_label']]
df_meta = df_meta.loc[df_meta.cluster_label=='glia', :]

# glia cluster_ID list
cluster_IDs = np.unique(df_meta['cluster_ID'].values)

# get top cluster-specific genes
marker_genes_dir = './preprocessed/marker_genes_tsv'
marker_genes_list = []
dict_df_mch = {}
for cluster_ID in cluster_IDs:
	df = pd.read_table(
			os.path.join(marker_genes_dir, cluster_ID+'.tsv'), 
			header=0, 
			index_col=0)
	dict_df_mch[cluster_ID] = df
	marker_genes = {'cluster_ID': cluster_ID, 
					'names': df.index.tolist()}
	marker_genes_list.append(marker_genes)

df_marker_genes = pd.DataFrame(marker_genes_list)
print(df_marker_genes.head())

# generate heatmap of gene expression for each cluster
# generate corrcoef matrix
df_corr = pd.DataFrame(columns=df_marker_genes['cluster_ID'], index=df_rna.columns)
df_pvalue = pd.DataFrame(columns=df_marker_genes['cluster_ID'], index=df_rna.columns)
for i, row in df_marker_genes.iterrows():
	# for each cluster, look at different genes

	# filter out top genes that are also in rna-seq list
	# top_genes = [gene if (gene in df_rna.index) for gene in row.names]
	top_genes = []
	for gene in row.names:
		if gene in df_rna.index:
			top_genes.append(gene)

	# get top 20 
	if len(top_genes) > num_cutoff:
		top_genes = top_genes[:num_cutoff]
	else:
		pass

	## select data from the df_rna

	# this doesn't keep the order of top genes
	# df_markers = df_rna.loc[df_rna.index.isin(top_genes), :] 
	# this keeps the order	
	idx = df_rna.index.get_indexer(top_genes)
	df_markers = df_rna.iloc[idx, :]

	# apply zscore to df_markers
	# df_markers = df_markers.apply(zscore, axis=1)
	
	#
	mch = dict_df_mch[row.cluster_ID].loc[top_genes, row.cluster_ID+'_mcc']	
	assert mch.index.tolist() == df_markers.index.tolist()
	print(mch.shape)
	print(df_markers.shape)
	corrs_res = [spearmanr(mch, df_markers[col].values)
				for col in df_rna.columns]	
	corrs = [item[0] for item in corrs_res]	
	pvalues = [item[1] for item in corrs_res]	

	df_corr[row.cluster_ID] = corrs
	# zscore for each column
	df_corr = df_corr.apply(zscore, axis=0)
	df_pvalue[row.cluster_ID] = pvalues 


	# generate heat map
	# fig, ax = plt.subplots(figsize=(8,6))
	# sns.heatmap(df_markers, ax=ax, cmap='viridis',
	# 			cbar_kws={'label':'$\mathregular{zscore of log_{10}(FPKM) on rows}$'})
	# ax.set_xticklabels(ax.get_xticklabels(), rotation=15)
	# ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
	# ax.set_title('Gene expressions of %s specific hypo-mCH genes' % row.cluster_ID)
	# fig.tight_layout()
	# fig.savefig('./preprocessed/marker_genes_rna/%s.pdf' % row.cluster_ID)
	# print('Saved figure to ./preprocessed/marker_genes_rna/%s.pdf' % row.cluster_ID)
	# plt.show()



# generate heat map for df_coor
df_pvalue[df_pvalue>0.05] = np.nan 
df_pvalue[df_pvalue<0.05] = '*'
df_pvalue = df_pvalue.fillna('')

fig, ax = plt.subplots(figsize=(8,6))
sns.heatmap(df_corr, ax=ax, cmap='viridis',
			annot = df_pvalue, fmt='',
			cbar_kws={'label':'Zscore (across rows) of Spearmanr correlation'})
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
ax.set_title('Correlation between genebody normalized mCH and gene expression '
	+ '($\mathregular{\log_{10}(FPKM)}$) \n (top %d cluster specific genes)' % num_cutoff)

fig.tight_layout()
fig.savefig('./preprocessed/marker_genes_rna/corr_matrix.pdf')
print('./preprocessed/marker_genes_rna/corr_matrix.pdf')
plt.show()


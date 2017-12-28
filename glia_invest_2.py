#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from scipy.stats import zscore
from scipy.stats import spearmanr
from scipy.stats import entropy

from snmcseq import plot_utils

context = 'CG'

NORMAL_SELECTION = False 
PAIR_WISE_SELECTION = False # odd, might corrupted 
ENTROPY_SELECTION_HIGH = False 
ENTROPY_SELECTION_LOW = False 


PLOT_GENES = False 
PLOT_VALUES = False 
PLOT_LABELS = False 

HARD_CODED_GENES =False 
if HARD_CODED_GENES:
	gs = os.listdir('./preprocessed/2017-10-31-marker_genes_mCG_glia_selected')
	gs = [g.split('_')[2].replace('.pdf', '') for g in gs] 
	HARD_CODED_GENES = sorted(list(set(gs))) 
	print(len(HARD_CODED_GENES))
	print(HARD_CODED_GENES)

	# HARD_CODED_GENES = ['SLC1A2',
	# 					'GPR124',] 

DO_CORR = False 

FANCY_HEATMAP = False 
PLAIN_HEATMAP = False 
GLIA_ONLY = False 

NEW_BACKSPIN = True

# file names
input_mcc = './preprocessed/glia/human_hv1_hv2_gene_cluster_vCG_mcg_matrix.tsv'
input_meta = './metadata/human_hv1_hv2_metadata_cells.tsv' 
input_rna = './metadata/barreslab_rnaseq.tsv'
input_gene_id_name = './metadata/gene_id_to_names.tsv'
tSNE_fname = './preprocessed/human_hv1_hv2_CG_genebody_nmcc_tsne.tsv'
gene_mcc_dir = ('/cndd/Public_Datasets/single_cell_methylome/gene_level/human/'
		+ 'genebody_mCG_combined-hs1-hs2_by_gene')
input_cluster = './preprocessed/human_hv1_hv2_CG_genebody_nmcc_backspin.tsv'

OPC = 'Oligodendrocyte Precursor Cell'
NFO = 'Newly Formed Oligodendrocyte'
MO = 'Myelinating Oligodendrocytes'
EC = 'Endothelial Cells'

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

def set_cell_type_number(cluster_label):
	"""
	"""
	cell_types = ['exci', 'inhi', 'glia']
	for i, cell_type in enumerate(cell_types):
		if cluster_label == cell_type:
			return i 
	raise ValueError('Cell type not found!')


def	tsne_and_boxplots(df, tx='tsne_x', ty='tsne_y', tc='mCH', bx='cluster_ID', by='mCH',
					output=None, show=True, close=False, 
					b_ylim=None, title=None):
	"""
	boxplot and tSNE plot

	xlim, ylim is set to facilitate displaying glial clusters only

	"""
	fig, axs = plt.subplots(2,1,figsize=(6,8))

	ax = axs[0]
	im = ax.scatter(df[tx], df[ty], s=2, 
		c=plot_utils.mcc_percentile_norm(df[tc].values))
	# ax.set_xlim([-40, 40])
	# ax.set_ylim([40, 100])
	if title:
		ax.set_title(title)
	ax.set_xlabel('tsne_x')
	ax.set_ylabel('tsne_y')
	ax.set_aspect('equal')
	clb = plt.colorbar(im, ax=ax)
	clb.set_label(tc, rotation=270, labelpad=10)

	ax = axs[1]
	sns.boxplot(x=bx, y=by, data=df, ax=ax)
	ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
	if b_ylim:
		ax.set_ylim(b_ylim)

	fig.tight_layout()
	if output:
		# output_dir = './preprocessed/marker_genes_%s' % method
		# if not os.path.exists(output_dir):
		# 	os.makedirs(output_dir)

		# output_fname = os.path.join(output_dir, '%s_%s.pdf' % (cluster_ID, gene_name))
		fig.savefig(output)
		print('Saved to ' + output) 
	if show:
		plt.show()
	if close:
		plt.close(fig)

### --- load data
# gene-expression fpkm
df_rna = pd.read_table(input_rna, header=0, index_col='Gene symbol')
# apply upper case to all entries of gene symbols
df_rna.index = df_rna.index.str.upper() 
# cell types
cell_types = df_rna.columns.tolist()
cell_types.remove('Description')
df_rna = df_rna[cell_types]
cell_types = rename_celltypes(cell_types) 
df_rna.columns = cell_types 
# take 1+log10
df_rna += 1
df_rna = df_rna.applymap(np.log10)
## normalize across all types 
# df_rna = df_rna.divide(df_rna.mean(axis=1), axis=0)
# df_rna = df_rna.subtract(df_rna.mean(axis=1), axis=0)
# df_rna = df_rna.divide(df_rna.std(axis=1), axis=0)
# df_rna = df_rna.apply(zscore, axis=1)
# ----
# print(cell_types)
# print(df_rna.head())

# glia_clusters from metadata
df_meta = pd.read_table(input_meta, header=0)
df_meta = df_meta[['Sample', 'cluster_ID', 'cluster_label', 'mCG/CG']]
if NEW_BACKSPIN:
	df_meta = df_meta.drop('cluster_ID', axis=1)
	df_cluster = pd.read_table(input_cluster, header=0)
	df_meta = pd.merge(df_meta, df_cluster, left_on='Sample', right_on='Sample', how='inner')
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

# gene-cluster mCH
df_mcc = pd.read_table(input_mcc, header=0, index_col='id')
mean_mccs = df_mcc.mean_mcc
df_mcc.drop('mean_mcc', axis=1, inplace=True)
## normalize across all cells
# df_mcc = df_mcc.divide(mean_mccs, axis=0)
# df_mcc = df_mcc.subtract(mean_mccs, axis=0)
# df_mcc = df_mcc.divide(df_mcc.std(axis=1), axis=0)
# df_mcc = df_mcc.apply(zscore, axis=1)
# ----

# gene id to name
df_gene_id_name = pd.read_table(input_gene_id_name, header=0, index_col='geneID')
df_gene_id_name = df_gene_id_name[['geneName']]

### end of data loading -------



### get feature gene ids
# add geneName to df_mcc
# df_mcc = pd.merge(df_gene_id_name, df_mcc, how='inner', left_index=True, right_index=True)

# change df_rna index to gene id
df_rna = pd.merge(df_gene_id_name, df_rna, how='inner', left_on='geneName', right_index=True)
df_rna.drop('geneName', axis=1, inplace=True)


### select feature genes and store it in the dict 
dict_feature_genes_all = dict()

### normal selection method --
if NORMAL_SELECTION: 
	method = 'normal'
	num_cutoff = 15 
	dict_feature_genes = dict()
	feature_gene_ids = []

	df_mcc_glia = df_mcc[[item+'_mcc' for item in cluster_IDs_glia]]
	lowp = df_mcc_glia.quantile(q=0.1, axis=1)
	highp = df_mcc_glia.quantile(q=0.5, axis=1)

	for cluster_ID_i in cluster_IDs_glia:
		diffs = (highp - 
				df_mcc_glia.loc[df_mcc_glia[cluster_ID_i+'_mcc']<lowp, cluster_ID_i+'_mcc'])	
		diffs.sort_values(ascending=False, inplace=True)	
		cluster_spec_gene_ids = diffs.head(num_cutoff).index
		dict_feature_genes[cluster_ID_i] = df_gene_id_name.loc[cluster_spec_gene_ids, 'geneName']
	
	# store to all
	for cluster_ID, series in dict_feature_genes.items():
		if cluster_ID not in dict_feature_genes_all:
			dict_feature_genes_all[cluster_ID] = dict_feature_genes[cluster_ID]	
		else:
			dict_feature_genes_all[cluster_ID] = (dict_feature_genes_all[cluster_ID].append(
				dict_feature_genes[cluster_ID]))
## end of normal selection method --


## pair-wise selection method --
"""
Pair-wise selection selects out genes that are low < 2% and high > 80%  
across ALL 82 clusters. Concatenating them as features is a bit strange. 

Could take only data from glial clusters and try merging glial clusters only.
"""
if PAIR_WISE_SELECTION: 
	method = 'pair-wise'
	num_cutoff = 5 
	dict_feature_genes = dict()
	feature_gene_ids = []
	df_mcc_pct = df_mcc.rank(axis=1, pct=True)
	for cluster_ID_i in cluster_IDs_glia:
		# get mCH < 0.02 percentile	
		df_spec_i = df_mcc.loc[df_mcc_pct[cluster_ID_i+'_mcc']<0.02, :]	
		cluster_spec_gene_ids = []
		for cluster_ID_j in cluster_IDs_glia:
			if cluster_ID_j != cluster_ID_i:
				# get mCH > 0.8 percentile
				df_spec_ij = df_spec_i.loc[df_mcc_pct[cluster_ID_j+'_mcc']>0.8, 
											[cluster_ID_i+'_mcc', cluster_ID_j+'_mcc']]
				# sort according to diff	
				df_spec_ij['diff'] = df_spec_ij[cluster_ID_j+'_mcc'] - df_spec_ij[cluster_ID_i+'_mcc']
				df_spec_ij.sort_values(['diff'], ascending=False, inplace=True)
				# get top marker gene ids
				spec_gene_ids = df_spec_ij.index.tolist()
				if len(spec_gene_ids) > num_cutoff:
					spec_gene_ids = spec_gene_ids[:num_cutoff]

				# print(cluster_ID_i, cluster_ID_j)
				# print(df_gene_id_name.loc[spec_gene_ids, 'geneName'])
				cluster_spec_gene_ids += spec_gene_ids
				feature_gene_ids += spec_gene_ids
		# 
		cluster_spec_gene_ids = list(set(cluster_spec_gene_ids)) 
		dict_feature_genes[cluster_ID_i] = df_gene_id_name.loc[cluster_spec_gene_ids, 'geneName']

	num_genes = len(feature_gene_ids)
	print(num_genes)
	feature_gene_ids = list(set(feature_gene_ids))
	df_mcc = df_mcc.loc[feature_gene_ids, :].sort_index()
	df_rna = df_rna.loc[feature_gene_ids, :].sort_index()
	assert df_mcc.index.tolist() == df_rna.index.tolist()
	num_genes = len(feature_gene_ids)
	print(num_genes)
## end of pair-wise selection method --

## entropy selection HIGH method --
if ENTROPY_SELECTION_HIGH:
	method = 'entropy_high'
	num_cutoff = 5 
	dict_feature_genes = dict()
	# glial cluster probability
	df_probs = df_mcc[[item+'_mcc' for item in cluster_IDs_glia]].divide(
				df_mcc[[item+'_mcc' for item in cluster_IDs_glia]].sum(axis=1), axis=0)

	entrops = df_probs.apply(entropy, axis=1)
	entrops.sort_values(inplace=True)
	dict_feature_genes['cluster_all'] = df_gene_id_name.loc[entrops.head(num_cutoff).index, 'geneName']

	# cluster-wise specific genes
	for cluster_ID in cluster_IDs_glia:
		Q_index = entrops -	df_probs[cluster_ID+'_mcc'].apply(np.log)	
		# top small indexes 
		Q_index.sort_values(inplace=True)
		# print(Q_index.head())
		# print(df_gene_id_name.loc[Q_index.head().index, 'geneName'])
		dict_feature_genes[cluster_ID] = df_gene_id_name.loc[Q_index.head().index, 'geneName']


	# store to all
	for cluster_ID, series in dict_feature_genes.items():
		if cluster_ID not in dict_feature_genes_all:
			dict_feature_genes_all[cluster_ID] = dict_feature_genes[cluster_ID]	
		else:
			dict_feature_genes_all[cluster_ID] = (dict_feature_genes_all[cluster_ID].append(
				dict_feature_genes[cluster_ID]))
## end of entropy selection method --

## entropy selection LOW method --
if ENTROPY_SELECTION_LOW:
	method = 'entropy_low'
	num_cutoff = 20 
	dict_feature_genes = dict()
	# glial clusters 
	df_mcc_glia = df_mcc[[item+'_mcc' for item in cluster_IDs_glia]]
	# subtract
	df_probs = - df_mcc_glia.subtract(df_mcc_glia.max(axis=1)+0.001, axis=0)
	# normalize
	df_probs = df_probs.divide(df_probs.sum(axis=1), axis=0)

	entrops = df_probs.apply(entropy, axis=1)
	entrops.sort_values(inplace=True)
	dict_feature_genes['cluster_all'] = df_gene_id_name.loc[entrops.head(10).index, 'geneName']
	# top entropy as correlation matrix features!

	# cluster-wise specific genes
	for cluster_ID in cluster_IDs_glia:
		Q_index = entrops -	df_probs[cluster_ID+'_mcc'].apply(np.log)	
		# top small indexes
		Q_index.sort_values(inplace=True)
		# print(Q_index.head())
		# print(df_gene_id_name.loc[Q_index.head().index, 'geneName'])
		dict_feature_genes[cluster_ID] = df_gene_id_name.loc[Q_index.head(num_cutoff).index, 'geneName']

	# store to all
	for cluster_ID, series in dict_feature_genes.items():
		if cluster_ID not in dict_feature_genes_all:
			dict_feature_genes_all[cluster_ID] = dict_feature_genes[cluster_ID]	
		else:
			dict_feature_genes_all[cluster_ID] = (dict_feature_genes_all[cluster_ID].append(
				dict_feature_genes[cluster_ID]))
## end of entropy selection method --

### report feature gene names
# feature_gene_names = df_gene_id_name.loc[feature_gene_ids, 'geneName']
# print(feature_gene_names)

### ---

### make boxplot and tSNE mCH plot for each gene on those clusters only

df_tSNE = pd.read_table(tSNE_fname, header=0)

if PLOT_GENES:

	for cluster_ID, series in dict_feature_genes.items():
		for gene_id, gene_name in series.iteritems():
			gene_fname = os.path.join(gene_mcc_dir, gene_id+'_m%s.txt' % context)
			df_gene = pd.read_table(gene_fname, header=None, index_col=0)
			df_gene.columns=['cell', 'm%s' % context, 'normalized_m%s' % context]
			df_gene = pd.merge(df_gene, df_meta[['Sample', 'cluster_ID', 'cluster_label']], 
				how='right', left_on='cell', right_on='Sample')
			df_gene = pd.merge(df_tSNE, df_gene,
				how='right', left_on='cells', right_on='cell')
			df_gene = df_gene[df_gene.cluster_label=='glia']
			df_gene.sort_values('cluster_ID', inplace=True)

			# boxplot and tSNE plot
			# tsne_and_boxplots(data, tx='tsne_x', ty='tsne_y', tc='m%s' % context, 
								# bx='cluster_ID', by='m%s' % context)
			fig, axs = plt.subplots(2,1,figsize=(6,8))

			ax = axs[0]
			im = ax.scatter(df_gene.tsne_x, df_gene.tsne_y, s=2, 
				c=plot_utils.mcc_percentile_norm(df_gene['m%s' % context].values))
			# ax.set_xlim([-40, 40])
			# ax.set_ylim([40, 100])
			ax.set_title('Gene level m%s: %s (%s marker)' % (context, gene_name, cluster_ID))
			ax.set_xlabel('tsne_x')
			ax.set_ylabel('tsne_y')
			ax.set_aspect('equal')
			clb = plt.colorbar(im, ax=ax)
			clb.set_label('m%s level' % context, rotation=270, labelpad=10)

			ax = axs[1]
			sns.boxplot(x='cluster_ID', y='m%s' % context, data=df_gene, ax=ax)
			ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
			# ax.set_ylim([-0.01, 0.15])

			fig.tight_layout()
			output_dir = './preprocessed/marker_genes_%s' % method
			if not os.path.exists(output_dir):
				os.makedirs(output_dir)

			output_fname = os.path.join(output_dir, '%s_%s.pdf' % (cluster_ID, gene_name))
			fig.savefig(output_fname)
			print('Saved to ' + output_fname) 
			plt.close(fig)


### plot several genes
if HARD_CODED_GENES:
	for gene_name in HARD_CODED_GENES:
		gene_id = df_gene_id_name[df_gene_id_name.geneName==gene_name].index[0]
		gene_fname = os.path.join(gene_mcc_dir, gene_id+'_m%s.txt' % context)
		df_gene = pd.read_table(gene_fname, header=None, index_col=0)
		df_gene.columns=['cell', 'm%s' % context, 'normalized_m%s' % context]
		df_gene = pd.merge(df_gene, df_meta[['Sample', 'cluster_ID', 'cluster_label']], 
			how='inner', left_on='cell', right_on='Sample')
		df_gene = pd.merge(df_tSNE, df_gene,
			how='inner', left_on='cells', right_on='cell')
		df_gene = df_gene[df_gene.cluster_label=='glia']
		df_gene.sort_values('cluster_ID', inplace=True)
		# boxplot and tSNE plot
		tsne_and_boxplots(df_gene, tx='tsne_x', ty='tsne_y', tc='m%s' % context, 
							bx='cluster_ID', by='m%s' % context, 
							title='Gene level m%s: %s' % (context, gene_name),
							show=False, close=True,
							output='./preprocessed/marker_genes_hardcoded/%s.pdf' % (gene_name))

### ---



### plot global mCG
if PLOT_VALUES:

	df_global_mCG = pd.merge(df_tSNE, df_meta[['Sample', 'cluster_ID', 'cluster_label', 'mCG/CG']], 
		how='right', left_on='cells', right_on='Sample')
	df_global_mCG = df_global_mCG[df_global_mCG.cluster_label=='glia']
	df_global_mCG = df_global_mCG.sort_values('cluster_ID')
	tsne_and_boxplots(df_global_mCG, tx='tsne_x', ty='tsne_y', tc='mCG/CG', bx='cluster_ID', by='mCG/CG',
					output='./preprocessed/tsne_and_boxplots_global_mCG.pdf', show=True, close=False,
					b_ylim=None, title='Global mCG level (glial clusters)')

### ---

### generate corrcoef matrix
if DO_CORR:

	# intersect 2 gene list
	print(len(df_mcc.index))
	print(len(df_rna.index))
	feature_gene_ids = list(set(df_mcc.index) & set(df_rna.index))
	print(len(feature_gene_ids))
	if dict_feature_genes_all:
		feature_gene_ids_v2 = []	
		for cluster_ID, series in dict_feature_genes_all.items():
			feature_gene_ids_v2 += series.index.tolist()
		feature_gene_ids = list(set(feature_gene_ids) & set(feature_gene_ids_v2))

	# make sure we have the same list of ordered genes for both matrices	
	df_mcc = df_mcc.loc[feature_gene_ids, :].sort_index()
	df_rna = df_rna.loc[feature_gene_ids, :].sort_index()
	assert df_mcc.index.tolist() == df_rna.index.tolist()
	num_genes = len(feature_gene_ids)
	print(num_genes)

	df_corr = pd.DataFrame(columns=cluster_IDs, index=cell_types)
	df_pvalue = pd.DataFrame(columns=cluster_IDs, index=cell_types)
	for cluster_ID in cluster_IDs:
		# for each cluster, look at different cell types 

		corrs_res = [spearmanr(df_mcc[cluster_ID+'_mcc'].values, df_rna[cell_type].values)
					for cell_type in cell_types]	
		corrs = [item[0] for item in corrs_res]	
		pvalues = [item[1] for item in corrs_res]	

		df_corr[cluster_ID] = corrs
		df_pvalue[cluster_ID] = pvalues 
### --- 


## generate heat map for df_coor
### --- fancy heatmap
if FANCY_HEATMAP:
	# df_pvalue[df_pvalue>0.05] = np.nan 
	# df_pvalue[df_pvalue<0.05] = ''
	# df_pvalue = df_pvalue.fillna('')

	df_annot = pd.DataFrame(columns=df_corr.columns, index=df_corr.index)
	for col, idx in df_corr.idxmin().iteritems():
		print(col, idx)
		df_annot.loc[idx, col] = '*'
	df_annot = df_annot.fillna('')

	colors=['C0', 'C1' ,'C2']
	row_num = [set_cell_type_number(
		df_meta.loc[df_meta.cluster_ID==cluster_ID, 'cluster_label'].values[0])
		for cluster_ID in cluster_IDs]

	fig = plt.figure(figsize=(12,8))
	ax1 = plt.subplot2grid((20,10), (0,0), colspan=10, rowspan=9)
	ax2 = plt.subplot2grid((20,10), (9,0), colspan=10, rowspan=9)
	ax3 = plt.subplot2grid((20,10), (19,0), colspan=10, rowspan=1)

	ax = ax1
	sns.heatmap(df_corr, 
				ax=ax, cmap='viridis',
				annot = df_annot, fmt='',
				# linewidth = 0.01,
				xticklabels = False,
				cbar_kws={'label':'Spearmanr correlation'})
	# ax.set_xticklabels([], rotation=90)
	ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
	ax.set_title(
		('Correlation between genebody normalized m%s and gene expression'
		+'($\mathregular{\log_{10}(FPKM)}$) \n (%d genes)') % (context, num_genes))

	ax = ax2
	# zscore for each column
	sns.heatmap(df_corr.apply(zscore, axis=0), 
				ax=ax, cmap='viridis',
				annot = df_annot, fmt='',
				# linewidth = 0.01,
				xticklabels = False,
				cbar_kws={'label':'Zscore (across rows)'})
	# ax.set_xticklabels([], rotation=90)
	ax.set_yticklabels(ax.get_yticklabels(), rotation=0)

	ax = ax3
	sns.heatmap(pd.DataFrame(row_num).transpose(), 
				ax=ax, 
				cmap=colors,
				# linewidth = 0.01,
				xticklabels = False,
				yticklabels = False,
				#annot = df_annot, fmt='',
				cbar_kws={'label': 'Exci, Inhi, Glia', 
						'ticks': [],
						}
				)
	# ax.set_xticklabels([], rotation=90)
	ax.set_yticklabels([], rotation=0)
	ax.set_xlabel('Clusters')

	# fig.tight_layout()
	output_fname = './preprocessed/corr_matrix_v2_%d.pdf' % num_genes
	fig.savefig(output_fname)
	print(output_fname)
	plt.show()
### ----


### --- plain heatmap
if PLAIN_HEATMAP:
	# df_value[df_pvalue>0.05] = np.nan 
	# df_pvalue[df_pvalue<0.05] = ''
	# df_pvalue = df_pvalue.fillna('')

	df_annot = pd.DataFrame(columns=df_corr.columns, index=df_corr.index)
	for col, idx in df_corr.idxmin().iteritems():
		print(col, idx)
		df_annot.loc[idx, col] = '*'
	df_annot = df_annot.fillna('')


	fig, axs = plt.subplots(2,1,figsize=(12,8))
	ax = axs[0]
	sns.heatmap(df_corr, 
				ax=ax, cmap='viridis',
				annot = df_annot, fmt='',
				# linewidth = 0.01,
				xticklabels = False,
				cbar_kws={'label':'Spearmanr correlation'})
	# ax.set_xticklabels([], rotation=90)
	ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
	ax.set_title(('Correlation between genebody normalized m%s and gene expression '
		+ '($\mathregular{\log_{10}(FPKM)}$) \n (%d genes)') % (context, num_genes))

	ax = axs[1]
	# zscore for each column
	sns.heatmap(df_corr.apply(zscore, axis=0), 
				ax=ax, cmap='viridis',
				annot = df_annot, fmt='',
				# linewidth = 0.01,
				# xticklabels = False,
				cbar_kws={'label':'Zscore (across rows)'})
	ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
	ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
	# ax.set_title('Correlation between genebody normalized mCH and gene expression '
	# 	+ '($\mathregular{\log_{10}(FPKM)}$) \n (All genes)')
	### ----


	# fig.tight_layout()
	output_fname = './preprocessed/corr_matrix_v2_glia_only_%d.pdf' % num_genes
	fig.savefig(output_fname)
	print(output_fname)
	plt.show()


## get cluster cell numbers

df_cluster_count = df_meta.groupby('cluster_ID').count()
df_cluster_count['cluster_ID'] = [int(item.strip('cluster_')) 
								for item in df_cluster_count.index.values]
df_cluster_count = df_cluster_count[['cluster_ID', 'Sample']].sort_values('cluster_ID')
df_cluster_count.columns = ['cluster_ID', 'Number of cells']

fig, ax = plt.subplots(figsize=(12,6))
colors = ['C0', 'C1', 'C2']
row_colors = [colors[set_cell_type_number(
	df_meta.loc[df_meta.cluster_ID==cluster_ID, 'cluster_label'].values[0])]
	for cluster_ID in cluster_IDs]
# bar_colors
sns.barplot(x='cluster_ID', y='Number of cells', data=df_cluster_count,
			ax=ax,
			# color=['C0', 'C1'],
			#palette=row_colors
			)	
ax.set_title("Number of cells in all clusters")
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
fig.savefig('./preprocessed/barplot_cluster_cell_numbers.pdf')
plt.show()
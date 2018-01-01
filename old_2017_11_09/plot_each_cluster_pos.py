#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from snmcseq.utils import set_cell_type_number
# import snmcseq.plot_utils

NEW_BACKSPIN = True

### load data
# tsne_coordinates
input_tsne = './preprocessed/human_hv1_hv2_CH_genebody_nmcc_tsne.tsv'
df_tsne = pd.read_table(input_tsne, header=0)

# metadata
input_meta = './metadata/human_hv1_hv2_metadata_cells.tsv' 
df_meta = pd.read_table(input_meta, header=0)
df_meta = df_meta[['Sample', 'cluster_ID', 'cluster_label', 'Batch']]
if NEW_BACKSPIN:
	input_cluster = './preprocessed/human_hv1_hv2_CG_genebody_nmcc_backspin.tsv'
	df_meta = df_meta.drop('cluster_ID', axis=1)
	df_cluster = pd.read_table(input_cluster, header=0)
	df_meta = pd.merge(df_meta, df_cluster, left_on='Sample', right_on='Sample', how='inner')
cluster_IDs = np.unique(df_meta.cluster_ID.values)


df_info = pd.merge(df_tsne, df_meta, how='inner', left_on='cells', right_on='Sample')

df_info_glia = df_info[df_info.cluster_label=='glia']
cluster_IDs_glia = np.unique(df_meta[df_meta.cluster_label=='glia'].values)

colors = ['C0', 'C1', 'C3']
for cluster_ID in cluster_IDs_glia:
	cell_type = df_meta.loc[df_meta.cluster_ID==cluster_ID, 'cluster_label'].values[0]

	color_col = [colors[set_cell_type_number(cell_type)] if item==cluster_ID else 'gray' 
			for item in df_info_glia.cluster_ID]	

	fig, ax = plt.subplots(figsize=(8,8))
	ax.scatter(df_info_glia.tsne_x, df_info_glia.tsne_y, 
			c=color_col, 
			s=8)
	# ax.set_xlim([-40, 40])
	# ax.set_ylim([40, 100])
	ax.set_aspect('equal')
	ax.set_xlabel('tSNE_x')
	ax.set_ylabel('tSNE_y')
	ax.set_title('%s' % cluster_ID)
	fig.savefig('./preprocessed/cluster_loc/%s.pdf' % cluster_ID)
	print('Saved to ./preprocessed/cluster_loc/%s.pdf' % cluster_ID)
	plt.show()
	# plt.clf()

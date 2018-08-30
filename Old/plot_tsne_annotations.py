#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt

from snmcseq_utils import plot_tsne_labels

celltype_res = pd.read_table('./data/cluster/cluster_MB_v1_MB_EA_MB_EB/cluster_MB_v1_MB_EA1_annotation.tsv')
# tSNE plot as confirmation
df_tsne = pd.read_table('./data/cluster/cluster_MB_v1_MB_EA1/cluster_MB_v1_MB_EA1_tsne.tsv')

# plot
df_plot = pd.merge(df_tsne, celltype_res, left_on='cluster_ID', right_on='cluster_ID', how='left')
plot_tsne_labels(df_plot, tc='cluster_annotation', 
	title='tSNE of human MB_v1 and MB_EA1 samples', 
	figsize=(8,8), legend_mode=1,
	output='./results/tsne_clusters/tsne_cluster_MB_v1_MB_EA1_annotation.pdf'
	)

celltype_res = pd.read_table('./data/cluster/cluster_MB_v1/cluster_MB_v1.tsv')
celltype_res = celltype_res[['sample', 'cell_type']] 
# plot
df_plot = pd.merge(df_tsne, celltype_res, left_on='sample', right_on='sample', how='left')
plot_tsne_labels(df_plot, tc='cell_type', 
	title='tSNE of human MB_v1 and MB_EA1 samples', 
	figsize=(8,8), legend_mode=1,
	output='./results/tsne_clusters/tsne_cluster_MB_v1_cell_type.pdf'
	)

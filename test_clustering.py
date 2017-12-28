#!/usr/bin/env python3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA
from sklearn.decomposition import FastICA 
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans
from snmcseq_utils import plot_tsne_labels

# load data
input_file = './data/binc/binc_mCH_human_combined_100000_summary_normalized.tsv'
tsne_file = './data/tsne/tsne_perp50_binc_mCHmCG_human_combined_100000_summary_normalized.tsv'
meta_file = './data/metadata/metadata_human_combined_updated.tsv' 
cluster_file = './data/backspin/clusters_binc_mCHmCG_human_combined_100000_summary_normalized.tsv'

df = pd.read_table(input_file)
print(df.shape)
print(df.columns.values[:5])
samples = [col[:-len('_mcc')] for col in df.columns]

df_tsne = pd.read_table(tsne_file, index_col='cells')

df_meta = pd.read_table(meta_file, index_col='Sample')
df_cluster = pd.read_table(cluster_file, index_col='Sample')
print(df_meta.shape)

### ----

# # clustering methods

# # k-means
# res = KMeans(n_clusters=60).fit_predict(df.T)

# print(res)
# print(res.shape)
# print(type(res))
# df_cluster = pd.DataFrame()
# df_cluster['Sample'] = samples 
# df_cluster = df_cluster.set_index('Sample')
# df_cluster['cluster_ID'] = res
# df_plot = pd.merge(df_tsne, df_cluster, left_index=True, right_index=True)



# backSPIN
df_plot = pd.merge(df_tsne, df_cluster, left_index=True, right_index=True)

# plot tsne clustering results
# plot_tsne_labels(df_plot, tc='cluster_ID', 
# 	title='tSNE of human MB_v1, MB_EA, MB_EB samples', 
# 	figsize=(8,8), legend_mode=1,
# 	output='./results/tsne_clusters/tsne_cluster_MB_v1_MB_EA_MB_EB_clusters.pdf'
# 	)


# plot individual clusters
for cluster_ID in np.unique(df_plot.cluster_ID.values):

	color_col = ['C0' if item==cluster_ID else 'gray' 
			for item in df_plot.cluster_ID]	

	fig, ax = plt.subplots(figsize=(8,8))
	ax.scatter(df_plot.tsne_x, df_plot.tsne_y, 
			c=color_col, 
			s=8)
	# ax.set_xlim([-40, 40])
	# ax.set_ylim([40, 100])
	# ax.set_aspect('equal')
	ax.set_aspect('auto')
	ax.set_xlabel('tSNE_x')
	ax.set_ylabel('tSNE_y')
	ax.set_title('%s' % cluster_ID)
	fig.savefig('./results/tsne_clusters/tsne_cluster_MB_v1_MB_EA_MB_EB/%s.pdf' % cluster_ID)
	print('Saved to ./results/tsne_clusters/tsne_cluster_MB_v1_MB_EA_MB_EB/%s.pdf' % cluster_ID)
	# plt.show()
	plt.clf()



# plot pca
# pcs = PCA(n_components=5).fit_transform(df.T)
# df_pcs = pd.DataFrame(pcs, index=samples)
# df_plot = pd.merge(df_meta, df_pcs, left_index=True, right_index=True)
# fig, ax = plt.subplots()
# for label, df_sub in df_plot.groupby('Biosample'):
# 	ax.scatter(df_sub[0], df_sub[2], s=1, label=label)
# ax.legend()
# plt.show()


# # plot ica
# N = 10
# ics = FastICA(n_components=N).fit_transform(df.T)
# df_ics = pd.DataFrame(ics, index=samples)
# df_plot = pd.merge(df_tsne, df_ics, left_index=True, right_index=True)
# for i in range(N):
# 	fig, ax = plt.subplots()
# 	ax.scatter(df_plot.tsne_x, df_plot.tsne_y, s=1, c=df_plot[i])
# 	# ax.legend()
# 	plt.show()



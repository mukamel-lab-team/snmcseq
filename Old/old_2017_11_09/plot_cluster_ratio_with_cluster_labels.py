#!/usr/bin/python3
"""
Generate a bar plot for cluster ratio labeled by cluster label (exci, inhi, glia)

"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def label_to_int(label):
	if label == 'exci':
		return 0
	elif label == 'inhi':
		return 1
	elif label == 'glia':
		return 2
	else:
		return 'nan'

### setup 
output_fname = './preprocessed/barplot_human_hv1_hv2_cluster_ratios.pdf'
output_fname_2 = './preprocessed/barplot_human_hv1_hv2_exci_inhi_glia_ratios.pdf'
output_fname_3 = './preprocessed/boxplot_cell_propotion_in_clusters.pdf'

### load data 
cluster_fname = './preprocessed/human_hv1_hv2_CH_100000_clusters.tsv'
df = pd.read_table(cluster_fname, header=0)

cluster_labels_fname = './preprocessed/human_hv1_hv2_CH_100000_cluster_labels.tsv'
# format: cluster_label/cluster_num
df_cl = pd.read_table(cluster_labels_fname, header=0)

### get info
df['biosample'] = [int(item.split('_')[0].strip('hv_'))
						for item in df['sample'].tolist()]
df['cluster_ID'] = [int(item.strip('cluster_'))
						for item in df['cluster_ID'].tolist()]
clusters = sorted(np.unique(df['cluster_ID'].values).tolist())
biosamples = sorted(np.unique(df['biosample'].values).tolist())
n_clusters = len(clusters)
n_biosamples = len(biosamples)

# get label info



df_ref = pd.DataFrame()
df_ref['cluster_ID'] = clusters	
df_ref = pd.merge(df_ref, df_cl, how='left', 
	left_on='cluster_ID', right_on='cluster_num', sort=True)
df_ref = df_ref[['cluster_ID', 'cluster_label']]

df_plot = df_ref.copy()

for biosample, df_bio in df.groupby('biosample'):
	count = df_bio.groupby('cluster_ID').count()
	count['cluster_ID'] = count.index.values
	df_new = pd.merge(df_ref, count, how='left',
				left_on='cluster_ID', right_on='cluster_ID', sort=True)
	df_new = df_new.fillna(0)
	# format: 1/2 (int) biosample
	df_plot[biosample] = df_new['sample'].values

# normalize
total_counts = np.zeros(n_clusters)
for biosample in biosamples:
	total_counts += df_plot[biosample].values

for biosample in biosamples:
	df_plot[biosample] = df_plot[biosample]/total_counts
# sort
df_plot['total'] = total_counts
df_plot.sort_values([biosamples[0]], ascending=False, inplace=True)
# format: cluster_ID/cluster_label/1(biosample)/2(biosample)/total


colors = ['b', 'r', 'g']
patterns = ['', '//', '/', 'o']
df_plot['cluster_color'] = [colors[label_to_int(label)] 
							for label in df_plot.cluster_label] 
# print(df_plot)

# #cluster ratio plot
# fig, axs = plt.subplots(2,1, figsize=(10,6))	
# ax = axs[0]
# offset = 0.0
# index = offset + np.arange(n_clusters, dtype=float)
# bar_width = 1.0 
# heights = np.zeros(n_clusters)
# for i, biosample in enumerate(biosamples):
# 	ax.bar(index, df_plot[biosample].values, bar_width, heights, 
# 		color=df_plot['cluster_color'],
# 		hatch=patterns[i],
# 		label='human_'+str(biosample))	
# 	heights += df_plot[biosample].values 
# ax.set_xlabel('Cluster number')
# ax.set_ylabel('Cell number propotion')
# handles, labels = ax.get_legend_handles_labels()
# ax.set_xticks(index+0.5)
# xticklabels = [item if i%1==0 else '' 
# 			for i, item in enumerate(df_plot.cluster_ID.tolist())]
# ax.set_xticklabels(xticklabels, rotation=90, fontsize='x-small')
# ax.legend()

# ax = axs[1]
# ax.bar(index, df_plot['total'].values, bar_width,  
# 		color=df_plot['cluster_color'],
# 		label='total number')	
# ax.set_xlabel('Cluster number')
# ax.set_ylabel('Total number of cells')
# handles, labels = ax.get_legend_handles_labels()
# # hard coded labels
# labels = ['Excitatory', 'Inhibitory', 'Glia']
# handles = [handles[0][0], handles[0][40], handles[0][-1]]
# ax.legend(handles, labels)
# # ax.set_title('Total number of cells in each cluster')
# ax.set_xticks(index+0.5)
# xticklabels = [item if i%3==0 else '' 
# 			for i, item in enumerate(df_plot.cluster_ID.tolist())]
# ax.set_xticklabels(xticklabels, rotation=90)

# fig.tight_layout()
# fig.savefig(output_fname)
# print("save file %s" % output_fname)


# data for a simplified plot

for biosample in biosamples:
	df_plot['num_'+str(biosample)] = df_plot['total'].values*df_plot[biosample].values
df_plot_simplified = df_plot[['cluster_label', 'total', 'num_1', 'num_2']].groupby('cluster_label').sum()
df_plot_simplified = df_plot_simplified.reindex(['exci', 'inhi', 'glia'])
# print(df_plot_simplified)
n_bars = 3

# # cluster ratio plot (simplified)
# fig, axs = plt.subplots(2,1, figsize=(8,6))	
# ax = axs[0]
# offset = 0.0
# index = offset + np.arange(n_bars, dtype=float)
# bar_width = 1.0 
# heights = np.zeros(n_bars)
# for i, biosample in enumerate(biosamples):

# 	barheight = df_plot_simplified['num_'+str(biosample)].values/df_plot_simplified['total'].values 
# 	print(barheight)
# 	ax.bar(index, 
# 		df_plot_simplified['num_'+str(biosample)].values/df_plot_simplified['total'].values, 
# 		bar_width, heights, 
# 		color=colors[:n_bars],
# 		hatch=patterns[i],
# 		label='human_'+str(biosample))	
# 	heights += np.reshape(barheight, (3,)) 
# ax.set_xlabel('Cluster number')
# ax.set_ylabel('Cell number propotion')
# handles, labels = ax.get_legend_handles_labels()
# ax.set_xticks(index+0.5)
# labels = ['Excitatory', 'Inhibitory', 'Glia']
# xticklabels = labels 
# ax.set_xticklabels(xticklabels, rotation=0)
# # xticklabels = [item if i%1==0 else '' 
# # 			for i, item in enumerate(df_plot.cluster_ID.tolist())]
# # ax.set_xticklabels(xticklabels, rotation=90, fontsize='x-small')
# ax.legend()



# ax = axs[1]
# ax.bar(index, df_plot_simplified['total'].values, bar_width,  
# 		color=colors[:n_bars],
# 		label='total number')	
# ax.set_xlabel('Cluster number')
# ax.set_ylabel('Total number of cells')
# handles, labels = ax.get_legend_handles_labels()

# # hard coded labels
# labels = ['Excitatory', 'Inhibitory', 'Glia']
# handles = [handles[0][0], handles[0][1], handles[0][2]]
# ax.legend(handles, labels)
# ax.set_title('Total number of cells in each cell type')
# ax.set_xticks(index+0.5)
# xticklabels = labels 
# ax.set_xticklabels(xticklabels, rotation=0)


# fig.tight_layout()
# fig.savefig(output_fname_2)
# print("save file %s" % output_fname_2)


### boxplot ---
data = []
labels = []
color_fills = []
fig, ax = plt.subplots(figsize=(8,6))
# b, r, g
colors = ['#1f77b4', '#d62728', '#2ca02c']
for cell_type, color in zip(['exci', 'inhi', 'glia'], colors):
	for biosample in biosamples:
		if biosample == 1:
			biosample_verbal = 'MFG'
		elif biosample == 2:
			biosample_verbal = 'BA10'
		labels.append('%s_%s' %(biosample_verbal, cell_type))
		color_fills.append(color)
		data.append(df_plot.loc[df_plot['cluster_label']==cell_type, biosample]) 

# bplot = ax.boxplot(data, 
# 			labels=labels,
# 			# positions = [1,1.5],
# 			patch_artist=True, 
# 			# boxprops=dict(facecolor=colors[0], color='black'),
# 			boxprops=dict(color='black'),
# 			capprops=dict(color='black'),
# 			whiskerprops=dict(color='black'),
# 			flierprops=dict(color='black', markeredgecolor='black'),
# 			medianprops=dict(color='black'))

# for patch, color in zip(bplot['boxes'], color_fills):
# 	patch.set_facecolor(color)

for i, (data_col, color) in enumerate(zip(data, color_fills)): 
	ax.scatter(i+np.zeros(len(data_col)), data_col, color=color)

ax.set_xticklabels(['']+labels, rotation=30)
ax.set_ylabel('Cell number proportion')
ax.set_ylim([-0.1,1.1])
fig.tight_layout()

# sns.boxplot(x='cluster_label', y=1, data=df_plot)
# sns.boxplot(x='cluster_label', y=2, data=df_plot)

fig.savefig(output_fname_3)
print("save file %s" % output_fname_3)

plt.show()	

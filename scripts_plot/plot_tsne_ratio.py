#!/usr/bin/python3

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sn
import time 


def plot_tsne_ratio(cluster_fname, 
				output_fname=None, title=None, show=None):
	"""
	"""
	df = pd.read_table(cluster_fname, header=0) 
	df['biosample'] = [int(item.split('_')[0].strip('hv_'))
							for item in df['sample'].tolist()]
	df['cluster_ID'] = [int(item.strip('cluster_'))
							for item in df['cluster_ID'].tolist()]

	clusters = sorted(np.unique(df['cluster_ID'].values).tolist())
	biosamples = sorted(np.unique(df['biosample'].values).tolist())
	n_clusters = len(clusters)
	n_biosamples = len(biosamples)

	df_ref = pd.DataFrame()
	df_ref['cluster_ID'] = clusters	
	df_plot = pd.DataFrame()
	df_plot['cluster_ID'] = clusters	

	for biosample, df_bio in df.groupby('biosample'):
		count = df_bio.groupby('cluster_ID').count()
		count['cluster_ID'] = count.index.values
		df_new = pd.merge(df_ref, count, how='left',
					left_on='cluster_ID', right_on='cluster_ID', sort=True)
		df_plot[biosample] = df_new['sample'].values

	# # format: cluster_ID/1/2 (int) biosample
	# fig, ax = plt.subplots()	
	# index = np.arange(n_clusters, dtype=float)
	# bar_width = 0.5
	# colors = ['b', 'r', 'y']
	# for i, biosample in enumerate(biosamples):
	# 	ax.bar(index, df_plot[biosample], bar_width, 
	# 		color=colors[i%len(colors)], label='human_'+str(biosample))	
	# 	index += bar_width
	# ax.legend()

	# format: cluster_ID/1/2 (int) biosample
	fig, ax = plt.subplots()	
	index = np.arange(n_clusters, dtype=float)
	bar_width = 0.9
	colors = ['b', 'r', 'y', 'g']
	heights = np.zeros(n_clusters)
	for i, biosample in enumerate(biosamples):
		ax.bar(index, df_plot[biosample].values, bar_width, heights, 
			color=colors[i%len(colors)], label='human_'+str(biosample))	
		heights += df_plot[biosample].values 
	ax.set_xlabel('Cluster number')
	ax.set_ylabel('Cell count')
	ax.legend()

	if title:
		ax.set_title(title)
	if output_fname:
		fig.savefig(output_fname)
		print("save file %s" % output_fname)
	if show:	
		plt.show()	



def create_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--input', required=True, help="input tsne results")

	parser.add_argument('-o', '--output', help="output tsne results")
	parser.add_argument('-t', '--title', help='title')
	parser.add_argument('-s', '--show', help='show plot', action="store_true")
	return parser

if __name__ == '__main__':

	parser = create_parser()
	args = parser.parse_args() 
	ti = time.time()

	plot_tsne_ratio(args.input, output_fname=args.output, 
					title=args.title, show=args.show)

	tf = time.time()
	print('time: %s' % (tf-ti))

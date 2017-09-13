#!/usr/bin/env python3

# tsne scatter plot from tsne matrix

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import argparse
import glob
import time



def plot_tsne_label(input_fname, meta_fname, label_col, merge_col='Sample', 
				output_fname=None, title=None, show=False, legend_mode=0):
	"""
	plot tsne colored by label from meta_fname
	"""
	# format: tsne_x, tsne_y, cells
	df = pd.read_table(input_fname, header=0)
	df_meta = pd.read_table(meta_fname, header=0)
	# format: Sample, Label_col 
	df_meta = df_meta[[merge_col, label_col]]	
	df_merged = pd.merge(df, df_meta, left_on='cells', right_on=merge_col, how='left', sort=True)

	fig, ax = plt.subplots()
	for label, df_sub in df_merged.groupby(label_col):
		ax.plot(df_sub['tsne_x'], df_sub['tsne_y'], 'o', label=label, markersize=1)
		
	ax.set_xlabel('tsne_x')
	ax.set_ylabel('tsne_y')
	if legend_mode == 0:
		ax.legend()
	elif legend_mode == 1:
		# Shrink current axis's height by 10% on the bottom
		box = ax.get_position()
		ax.set_position([box.x0, box.y0 + box.height * 0.1,
		                 box.width, box.height * 0.9])
		ax.legend(loc='upper center', bbox_to_anchor=(0.5, 0.06),
	          ncol=5, fancybox=False, shadow=False)	
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
	parser.add_argument('-m', '--metadata', required=True, help="meta data")
	parser.add_argument('-l', '--label_col', required=True, help="column in meta data as the label")
	parser.add_argument('-o', '--output', help="output tsne results")
	parser.add_argument('-t', '--title', help='title')
	parser.add_argument('-s', '--show', help='show plot', action="store_true")
	parser.add_argument('-lm', '--legend_mode', type=int, default=0, help="column in meta data as the label")
	return parser


if __name__ == '__main__':

	parser = create_parser()
	args = parser.parse_args() 
	ti = time.time()
	
	# plot_tsne_label('../tsne/human_combined_MB_EB_CH_100000_out_v1_25.tsv', 
	# 	'../metadata/MB_EB_metadata_cells.tsv', 'Layer', title='Deep v.s. Superficial', show=True)

	# plot_tsne_label('../tsne/human_combined_hv1_hv2_CH_100000_out_v1_10.tsv', 
	# 	'../metadata/human_hv1_hv2_metadata_cells.tsv', 'Batch', title='tsne_2_human', show=True)

	plot_tsne_label(args.input, args.metadata, args.label_col, output_fname=args.output,
				 title=args.title, show=args.show, legend_mode=args.legend_mode) 
	
	tf = time.time()
	print('time: %s' % (tf-ti))
 
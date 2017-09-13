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


def set_value_by_percentile(this, lo, hi):
	"""set `this` below or above percentiles to given values
	this (float)
	lo(float)
	hi(float)
	"""
	if this < lo:
		return lo
	elif this > hi:
		return hi
	else:
		return this

def mcc_percentile_norm(mcc, low_p=5, hi_p=95):
	"""
	set values above and below specific percentiles to be at the value of percentiles 

	args: mcc, low_p, hi_p	

	return: normalized mcc levels
	"""
#	mcc_norm = [np.isnan(mcc) for mcc_i in list(mcc)]
	mcc_norm = np.copy(mcc)
	mcc_norm = mcc_norm[~np.isnan(mcc_norm)]

	lo = np.percentile(mcc_norm, low_p)
	hi = np.percentile(mcc_norm, hi_p)

	mcc_norm = [set_value_by_percentile(mcc_i, lo, hi) for mcc_i in list(mcc)]
	mcc_norm = np.array(mcc_norm)

	return mcc_norm


def plot_tsne_value(input_fname, meta_fname, value_col, 
				output_fname=None, title=None, show=False):
	"""
	plot tsne colored by label from meta_fname
	"""
	# format: tsne_x, tsne_y, cells
	df = pd.read_table(input_fname, header=0)
	df_meta = pd.read_table(meta_fname, header=0)
	# format: Sample, Layer
	df_meta = df_meta[['Sample', value_col]]	
	df_merged = pd.merge(df, df_meta, left_on='cells', right_on='Sample', how='left', sort=True)

	fig, ax = plt.subplots()
	colors = df_merged[value_col].values 
	colors = mcc_percentile_norm(colors)
	size = [1]*df_merged.shape[0]

	im = ax.scatter(df_merged['tsne_x'], df_merged['tsne_y'], c=colors, s=size)
	fig.colorbar(im)

	ax.set_xlabel('tsne_x')
	ax.set_ylabel('tsne_y')
	
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
	parser.add_argument('-l', '--value_col', required=True, help="column in meta data as the label")
	parser.add_argument('-o', '--output', help="output tsne results")
	parser.add_argument('-t', '--title', help='title')
	parser.add_argument('-s', '--show', help='show plot', action="store_true")
	return parser


if __name__ == '__main__':

	parser = create_parser()
	args = parser.parse_args() 
	ti = time.time()
	
	# plot_tsne_label('../tsne/human_combined_MB_EB_CH_100000_out_v1_25.tsv', 
	# 	'../metadata/MB_EB_metadata_cells.tsv', 'Layer', title='Deep v.s. Superficial', show=True)

	# plot_tsne_label('../tsne/human_combined_hv1_hv2_CH_100000_out_v1_10.tsv', 
	# 	'../metadata/human_hv1_hv2_metadata_cells.tsv', 'Batch', title='tsne_2_human', show=True)

	plot_tsne_value(args.input, args.metadata, args.value_col, output_fname=args.output,
				 title=args.title, show=args.show) 
	
	tf = time.time()
	print('time: %s' % (tf-ti))

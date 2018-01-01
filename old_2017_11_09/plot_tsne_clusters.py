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


def plot_tsne_label_v2(input_fname, ref_fname, label_col='cluster_ID', merge_col='Sample', 
				output_fname=None, title=None, show=False, legend_mode=0):
	"""
	plot tsne colored by label from df_ref
	df_ref format: [merge_col, label_col]
	"""

	# format: tsne_x, tsne_y, cells
	df = pd.read_table(input_fname, header=0)
	df_ref = pd.read_table(ref_fname, header=0)
	# format: sample, Label_col 
	df_ref = df_ref[[merge_col, label_col]]	
	df_merged = pd.merge(df, df_ref, left_on='cells', right_on=merge_col, how='left', sort=True)

	fig, ax = plt.subplots()
	fig.set_size_inches(11, h=8.5)
	for label, df_sub in df_merged.groupby(label_col):
		ax.plot(df_sub['tsne_x'], df_sub['tsne_y'], 'o', 
			label=label, markersize=4)
		
	if legend_mode == 0:
		ax.legend()
	elif legend_mode == 1:
		# Shrink current axis's height by 10% on the bottom
		box = ax.get_position()
		ax.set_position([box.x0, box.y0 + box.height * 0.1,
		                 box.width, box.height * 0.9])
		ax.legend(loc='upper center', bbox_to_anchor=(0.5, 0.06),
	          ncol=5, fancybox=False, shadow=False)	
		
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
	parser.add_argument('-r', '--refdata', required=True, help="reference data")
	parser.add_argument('-l', '--label_col', required=True, help="column in meta data as the label")
	parser.add_argument('--merge_col', required=True, help="column in ref data to merge the original ones")
	parser.add_argument('-o', '--output', help="output tsne results")
	parser.add_argument('-t', '--title', help='title')
	parser.add_argument('-s', '--show', help='show plot', action="store_true")
	parser.add_argument('-lm', '--legend_mode', type=int, default=0, help="column in meta data as the label")
	return parser


if __name__ == '__main__':

	parser = create_parser()
	args = parser.parse_args() 
	ti = time.time()
	
	plot_tsne_label_v2(args.input, args.refdata, args.label_col, 
				merge_col=args.merge_col, output_fname=args.output,
				title=args.title, show=args.show, 
				legend_mode=args.legend_mode) 
	
	tf = time.time()
	print('time: %s' % (tf-ti))
 
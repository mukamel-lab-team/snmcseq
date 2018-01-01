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


from snmcseq_utils import plot_tsne_labels

def plot_tsne_metalabel(input_fname, meta_fname, label_col, merge_col='Sample', 
				output_fname=None, title=None, show=False, close=False, legend_mode=0, figsize=(8,6),
				colors=['C0', 'C2', 'C1', 'C3', 'C4', 'C5', 'C6', 'C8', 'C9']):
	"""
	plot tsne colored by label from meta_fname
	"""
	# format: tsne_x, tsne_y, cells
	df = pd.read_table(input_fname)
	df_meta = pd.read_table(meta_fname)
	# format: Sample, Label_col 
	df_meta = df_meta[[merge_col, label_col]]	
	df_merged = pd.merge(df, df_meta, left_on='sample', right_on=merge_col, how='left', sort=True)

	plot_tsne_labels(df_merged, tx='tsne_x', ty='tsne_y', tc=label_col, 
                    legend_mode=0,
                    output=output_fname, show=show, close=close, 
                    t_xlim=None, t_ylim=None, title=title, figsize=figsize,
                    colors=colors)


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
	

	plot_tsne_metalabel(args.input, args.metadata, args.label_col, output_fname=args.output,
				 title=args.title, show=args.show, legend_mode=args.legend_mode) 
	
	tf = time.time()
	print('time: %s' % (tf-ti))
 

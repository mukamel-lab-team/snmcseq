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


from snmcseq_utils import plot_tsne_values

def plot_tsne_metavalues(input_fname, meta_fname, value_col, 
				output_fname=None, 
				s=2, title=None, show=False, close=False, figsize=(8,6)):
	"""
	plot tsne colored by label from meta_fname
	"""
	# format: tsne_x, tsne_y, cells
	df = pd.read_table(input_fname, header=0)
	df_meta = pd.read_table(meta_fname, header=0)
	# format: Sample, Layer
	df_meta = df_meta[['Sample', value_col]]	
	df_merged = pd.merge(df, df_meta, left_on='sample', right_on='Sample', how='left', sort=True)

	plot_tsne_values(df_merged, tx='tsne_x', ty='tsne_y', tc=value_col,
					s=s,
                    cbar_label=value_col,
                    output=output_fname, show=show, close=close, 
                    t_xlim=None, t_ylim=None, title=title, figsize=figsize)



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
	

	plot_tsne_metavalues(args.input, args.metadata, args.value_col, output_fname=args.output,
				 title=args.title, show=args.show) 
	
	tf = time.time()
	print('time: %s' % (tf-ti))

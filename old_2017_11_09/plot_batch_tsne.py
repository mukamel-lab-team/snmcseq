#!/usr/bin/env python3

# tsne scatter plot from tsne matrix

import numpy as np
import glob
import os
import matplotlib.pyplot as plt
import pandas as pd
import argparse
from collections import OrderedDict

def plot_batch_tsne(input_fnames, output_fname=None, titles=None):
	"""
	args: input_fnames (tsne output)

	output: a combined tsne plot with subplots 
	"""
	num = len(input_fnames)
	ncol = 3
	if num % ncol == 0:
		nrow = int(num/ncol)
	else:
		nrow = int(num/ncol) + 1 
	fig, axs = plt.subplots(nrow, ncol)
	for i, input_fname in enumerate(input_fnames):	
		df = pd.read_table(input_fname, header=0)
		# df = df[df['tsne_y']<100]
		ax = axs[int(i/3)][i%3]
		ax.plot(df['tsne_x'], df['tsne_y'], 'o', markersize=1)
		ax.set_xlabel('tsne_x')
		ax.set_ylabel('tsne_y')
		if titles:
			ax.set_title(titles[i])

	fig.tight_layout()
	if output_fname:
		fig.savefig(output_fname)
	plt.show()	

# def plot_batch_tsne_methyl(input_fnames, output_fname=None, titles=None):
# 	"""
# 	args: input_fnames (tsne output)

# 	output: a combined tsne plot with subplots 
# 	"""
# 	num = len(input_fnames)
# 	ncol = 3
# 	if num % ncol == 0:
# 		nrow = int(num/ncol)
# 	else:
# 		nrow = int(num/ncol) + 1 
# 	fig, axs = plt.subplots(nrow, ncol)
# 	for i, input_fname in enumerate(input_fnames):	
# 		df = pd.read_table(input_fname, header=0)
# 		# df = df[df['tsne_y']<100]
# 		ax = axs[int(i/3)][i%3]
# 		ax.plot(df['tsne_x'], df['tsne_y'], 'o', markersize=1)
# 		ax.set_xlabel('tsne_x')
# 		ax.set_ylabel('tsne_y')
# 		if titles:
# 			ax.set_title(titles[i])

# 	fig.tight_layout()
# 	if output_fname:
# 		fig.savefig(output_fname)
# 	plt.show()	
# def create_parser():
# 	parser = argparse.argumentparser()
# 	parser.add_argument('-i', '--input', required=true, help="input tsne results")
# 	parser.add_argument('-o', '--output', help="output tsne results")
# 	parser.add_argument('-t', '--title', help='title')
# 	return parser

if __name__ == '__main__':

	DIR = '/cndd/fangming/side_project/snmcseq_dev/tsne'
	input_fnames = glob.glob(
		os.path.join(DIR, 'human_combined_MB_EB_CH_100000_out_v1_*.tsv'))	
	titles = [os.path.basename(input_fname).split('.')[0] 
		for input_fname in input_fnames]
	perplexities = [int(title.split('_')[-1]) for title in titles]

	titles = ["perplexity = %d" % perplexity
		for perplexity in perplexities]
	dict_info = dict(zip(perplexities, zip(input_fnames, titles)))
	dict_info = OrderedDict(sorted(dict_info.items(), key=lambda t: t[0]))
	# print(dict_info)
	perplexities = list(dict_info.keys())
	items = list(dict_info.values())
	input_fnames = [item[0] for item in items]
	titles = [item[1] for item in items]
	# print(input_fnames)
	# print(titles)
	# print(perplexities)
	input_fnames = input_fnames[6:]
	titles = titles[6:]
	plot_batch_tsne(input_fnames, output_fname=None, titles=titles)

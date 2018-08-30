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

def plot_tsne_scatter(input_fname, output_fname=None, title=None):
	"""
	make basic tsne scatter plot
	"""
	df = pd.read_table(input_fname, header=0)

	fig, ax = plt.subplots()
	ax.plot(df['tsne_x'], df['tsne_y'], 'o', markersize=1)
	ax.set_xlabel('tsne_x')
	ax.set_ylabel('tsne_y')
	if title:
		ax.set_title(title)
	if output_fname:
		fig.savefig(output_fname)

	plt.show()	

# def get_sname_from_fname(fname):
# 	"""
# 	args: file name of a genebody methylation file

# 	return: sample name such as "Pool_2256_AD010_indexed_R1_mch_genebody.txt"
# 	"""
# 	sample_name = os.path.basename(fname).split('.')[0] 
# 	if sample_name.endswith('_mch_genebody'):
# 		sample_name = sample_name[:-len('_mch_genebody')]	

# 	return sample_name

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


def plot_tsne_scatter_gene_mCH(input_fname, gene_name, genebody_mch_dir, 
				output_fname=None, title=None, show=False):
	"""
	make a plot colored with marker gene methylation levels 
	"""
	df = pd.read_table(input_fname, header=0)

	# mcc = np.random.randn(df.shape[0]) 
	# mcc[-1] == 200

	# loop over all samples (cells) to get a list of mcc levels as color
	mcc = []
	for sample_name in df['cells']:	
		gene_mch_fname = glob.glob(os.path.join(genebody_mch_dir, sample_name+'*'))
		assert len(gene_mch_fname) == 1
		gene_mch_fname = gene_mch_fname[0]

		df_gene = pd.read_table(gene_mch_fname, header=0)
		df_gene = df_gene[df_gene['name'] == gene_name]
		mcc.append((df_gene['mc']/df_gene['c']).values[0])	

	size = [1]*df.shape[0]
	mcc_norm = mcc_percentile_norm(mcc)	

	fig, ax = plt.subplots()
	im = ax.scatter(df['tsne_x'], df['tsne_y'], c=mcc_norm, s=size)
	ax.set_xlabel('tsne_x')
	ax.set_ylabel('tsne_y')
	fig.colorbar(im)
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
	return parser

def create_parser_gene():
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--input', required=True, help="input tsne results")
	parser.add_argument('-g', '--gene', required=True, help="gene name")
	parser.add_argument('-gd', '--gene_dir', required=True, help="genebody mCH files")
	parser.add_argument('-o', '--output', help="output tsne results")
	parser.add_argument('-t', '--title', help='title')
	parser.add_argument('-s', '--show', help='show plot', action="store_true")
	return parser

if __name__ == '__main__':

	# input_fname = '/cndd/fangming/side_project/snmcseq_dev/tsne/tsne_out.tsv'
	# output_fname = '/cndd/fangming/side_project/snmcseq_dev/tsne/tsne_scatter_plot.pdf'
	# title = 'mCH_tsne_human_MB_EB'
	# parser = create_parser()
	# args = parser.parse_args() 
	# plot_tsne_scatter(args.input, output_fname=args.output, title=args.title)

	parser = create_parser_gene()
	args = parser.parse_args() 
	ti = time.time()
	plot_tsne_scatter_gene_mCH(args.input, args.gene, args.gene_dir, 
				output_fname=args.output, title=args.title,
				show=args.show)
	# plot_tsne_scatter_gene_mCH('./tsne/human_combined_MB_EB_CH_100000_out_v1_25.tsv', 
	# 	'GAD1', '/cndd/Public_Datasets/single_cell_methylome/gene_level/human/genebody_MB_EB')

	tf = time.time()
	print('time: %s' % (tf-ti))
 
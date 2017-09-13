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
	mcc_norm = np.copy(mcc)
	mcc_norm = mcc_norm[~np.isnan(mcc_norm)]

	lo = np.percentile(mcc_norm, low_p)
	hi = np.percentile(mcc_norm, hi_p)

	mcc_norm = [set_value_by_percentile(mcc_i, lo, hi) for mcc_i in list(mcc)]
	mcc_norm = np.array(mcc_norm)

	return mcc_norm

def gene_name_to_gene_id(gene_name, 
				gene_id_to_name_fname='/cndd/fangming/side_project/snmcseq_dev/metadata/gene_id_to_names.tsv'):
	"""
	"""	
	df_id_to_name = pd.read_table(gene_id_to_name_fname, header=0)
	gene_id = df_id_to_name[df_id_to_name['geneName']==gene_name]['geneID'].tolist()

	if len(gene_id) != 1:
		raise ValueError('Query error for gene ID')
	gene_id = gene_id[0]

	return gene_id

def plot_tsne_scatter_gene_mCH_v2(input_fname, gene_name, genebody_mch_dir, 
				gene_id_to_name_fname='/cndd/fangming/side_project/snmcseq_dev/metadata/gene_id_to_names.tsv', 
				output_fname=None, title=None, show=False, normalize=False, meta_fname='/cndd/fangming/side_project/snmcseq_dev/metadata/MB_EB_metadata_cells.tsv'):
	"""
	make a plot colored with marker gene methylation levels 
	"""
	# tsne df
	df = pd.read_table(input_fname, header=0)

	gene_id = gene_name_to_gene_id(gene_name, gene_id_to_name_fname=gene_id_to_name_fname)

	gene_fname = os.path.join(genebody_mch_dir, gene_id+'_mCH.txt')	
	if not os.path.isfile(gene_fname):
		raise ValueError('No such file: ' + gene_fname)
	# format: 0/1/2/3 -- gene_id, cell, mcc, mcc_norm
	df_gene = pd.read_table(gene_fname, header=None)
	df_merged = pd.merge(df, df_gene, 
		left_on='cells', right_on=1, sort=True)

	mcc = df_merged[2].tolist()

	if normalize:
		df_meta = pd.read_table(meta_fname, header=0)
		df_meta = df_meta[['Sample', 'mCH/CH']]	
		df_merged = pd.merge(df_merged, df_meta, 
					left_on='cells', right_on='Sample', sort=True)
		mcc = list(df_merged[2].values/df_merged['mCH/CH'].values)

	size = [1]*df.shape[0]
	mcc_norm = mcc_percentile_norm(mcc)	

	fig, ax = plt.subplots()
	im = ax.scatter(df_merged['tsne_x'], df_merged['tsne_y'], c=mcc_norm, s=size)
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

def plot_tsne_scatter_gene_mCH_v3(input_fname, gene_names, genebody_mch_dir, 
				gene_id_to_name_fname='/cndd/fangming/side_project/snmcseq_dev/metadata/gene_id_to_names.tsv', 
				output_fname=None, show=False):
	"""
	make a plot colored with marker gene methylation levels 
	"""
	# set up frames
	num = len(gene_names)
	ncol = 3
	if num % ncol == 0:
		nrow = int(num/ncol)
	else:
		nrow = int(num/ncol) + 1 

	fig, axs = plt.subplots(nrow, ncol)

	# tsne coords
	df = pd.read_table(input_fname, header=0)
	# loop over all genes
	for i, gene_name in enumerate(gene_names):
		# decide location
		ax = axs[int(i/ncol)][i%ncol]

		gene_id = gene_name_to_gene_id(gene_name, gene_id_to_name_fname=gene_id_to_name_fname)
		gene_fname = os.path.join(genebody_mch_dir, gene_id+'_mCH.txt')	
		if not os.path.isfile(gene_fname):
			raise ValueError('No such file: ' + gene_fname)
		# format: 0/1/2/3 -- gene_id, cell, mcc, mcc_norm
		df_gene = pd.read_table(gene_fname, header=None)
		df_merged = pd.merge(df, df_gene, 
			left_on='cells', right_on=1, sort=True)

		mcc = df_merged[2].tolist()

		size = [1]*df.shape[0]
		mcc_norm = mcc_percentile_norm(mcc)	

		im = ax.scatter(df_merged['tsne_x'], df_merged['tsne_y'], c=mcc_norm, s=size)
		ax.set_xlabel('tsne_x')
		ax.set_ylabel('tsne_y')
		fig.colorbar(im, ax=ax)
		ax.set_title(gene_name)
	fig.tight_layout()	

	if output_fname:
		fig.savefig(output_fname)
		print("save file %s" % output_fname)
	if show:	
		plt.show()	

def create_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--input', required=True, help="input tsne results")
	parser.add_argument('-g', '--gene', required=True, help="gene name")
	parser.add_argument('-gd', '--gene_dir', required=True, help="genebody mCH files")
	parser.add_argument('-o', '--output', help="output tsne results")
	parser.add_argument('-t', '--title', help='title')
	parser.add_argument('-n', '--normalize', help='normalize by global mCH level', action='store_true')
	parser.add_argument('-s', '--show', help='show plot', action="store_true")
	return parser


if __name__ == '__main__':

	parser = create_parser()
	args = parser.parse_args() 
	ti = time.time()
	# plot_tsne_scatter_gene_mCH_v2('../tsne/human_combined_MB_EB_CH_100000_out_v1_25.tsv', 
	# 	'GAD1', 
	# 	'/cndd/Public_Datasets/single_cell_methylome/gene_level/human/genebody_MB_EB_by_gene',
	# 	show=True)

	plot_tsne_scatter_gene_mCH_v2(args.input, 
								args.gene, 
								args.gene_dir,
								output_fname=args.output,
								title=args.title,
								normalize=args.normalize,
								show=args.show)

	tf = time.time()
	print('time: %s' % (tf-ti))
 

#!/usr/bin/env python3

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import argparse
import glob
import time

def plot_tsne_label_joint(input_fname,
				ref_fname, label_col, merge_col, 
				output_fname=None, title=None, show=False, legend_mode=0):
	"""
	plot tsne colored by label from meta_fname

	args:
		input_fname and refdata has to be hv1_hv2 tSNE coordinates and hv1_hv2 metadata
	"""
	# format: tsne_x, tsne_y, cells
	df = pd.read_table(input_fname, header=0)
	# format: merge_col2, label_col2	
	df_ref = pd.read_table(ref_fname, header=0)
	df_ref = df_ref[[merge_col, label_col]]

	df_merged = pd.merge(df, df_ref, left_on='cells', right_on=merge_col, how='left', sort=True)
	df_merged['biosample'] = [int(sample.split('_')[0].strip('hv')) 
							for sample in df_merged['cells'].tolist()]

	fig, ax = plt.subplots(figsize=(6,8))
	shapes = ['o', 's']
	# colors = ['C1', 'C2', 'C3']
	colors = ['b', 'g', 'r']
	for (i, (label, df_sub)) in enumerate(df_merged.groupby(label_col)):
		for biosample, df_sub_2 in df_sub.groupby('biosample'):
			# if biosample == 1:
			#	continue
			ax.plot(df_sub_2['tsne_x'], df_sub_2['tsne_y'],
				shapes[(biosample-1)%len(shapes)], 
				# color='C'+str(i%10), 
				color=colors[i%len(colors)], 
				label=label+'_human#'+str(biosample), markersize=5, markeredgecolor='black', alpha=0.9)

	ax.set_xlabel('tsne_x')
	ax.set_ylabel('tsne_y')
	ax.set_aspect('equal')
	if legend_mode == -1:
		pass
	elif legend_mode == 0:
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
	parser.add_argument('-r', '--refdata', required=True, help="ref data")
	parser.add_argument('-l', '--label_col', required=True, help="column in meta data as the label")
	parser.add_argument('--merge_col', required=True, help="sample column in ref data")

	parser.add_argument('-o', '--output', help="output tsne results")
	parser.add_argument('-t', '--title', help='title')
	parser.add_argument('-s', '--show', help='show plot', action="store_true")
	parser.add_argument('-lm', '--legend_mode', type=int, default=0, help="column in meta data as the label")
	return parser


if __name__ == '__main__':

	parser = create_parser()
	args = parser.parse_args() 
	ti = time.time()

	plot_tsne_label_joint(args.input,
				args.refdata, args.label_col, args.merge_col, 
				output_fname=args.output, title=args.title, show=args.show,
				legend_mode=args.legend_mode)

	tf = time.time()
	print('time: %s' % (tf-ti))
 

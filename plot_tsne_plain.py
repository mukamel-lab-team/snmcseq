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

def create_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--input', required=True, help="input tsne results")
	parser.add_argument('-o', '--output', help="output tsne results")
	parser.add_argument('-t', '--title', help='title')
	return parser


if __name__ == '__main__':
	
	parser = create_parser()
	args = parser.parse_args() 
	plot_tsne_scatter(args.input, output_fname=args.output, title=args.title)

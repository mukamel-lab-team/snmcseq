#!/usr/bin/env python3

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import ipdb
import numpy as np
import sys 
import argparse
import warnings
import matplotlib.pyplot as plt

# import mypy
# from collections import OrderedDict
# import random
# import math

GLIA_ONLY = False 
FLIP_AXES = True

def get_sample_names(df):
	"""
	"""
	col_names = df.columns.tolist()

	samples = []
	for col_name in col_names:
		if col_name.endswith('_c'):
			col_name = col_name[:-len('_c')] 
		elif col_name.endswith('_mc'):
			col_name = col_name[:-len('_mc')]
		else:
			raise 	
		samples.append(col_name)

	sample_names = list(set(samples))
	return sample_names


warnings.filterwarnings("ignore")

infile = './preprocessed/human_hv1_hv2_CH_genebody_nmcc.mat' 
outfile = './preprocessed/human_hv1_hv2_CH_genebody_nmcc_tsne_all.tsv'
# outfile_plot = './preprocessed/human_hv1_hv2_CG_genebody_nmcc_tsne.pdf'
perplexity = 25 
seed = 1

#######################
# LOADING DATA (preprocessed)
#######################
print('Loading data...')
df = pd.read_table(infile, header=0) 

# get glial cells
if GLIA_ONLY:
	df_meta = pd.read_table('./metadata/human_hv1_hv2_metadata_cells.tsv', 
		header=0, index_col='Sample')
	samples = df_meta[df_meta.cluster_label=='glia'].index.tolist()
	df = df[[sample+'_mcc' for sample in samples]]
print(df.shape)

#######################
# RUNNING PCA
#######################
print("Running PCA.")
pca = PCA(n_components=50)
sklearn_transf = pca.fit_transform(df.T)
# print(pca.explained_variance_ratio_)
sklearn_transf_PCA = sklearn_transf



#######################
# RUNNING TSNE
#######################
print("Running TSNE.")
num_components = 2
tsne = TSNE(n_components=num_components, init='pca', random_state=seed, perplexity=perplexity, verbose=3)
sklearn_transf = tsne.fit_transform(sklearn_transf_PCA)

print("Saving output to file.")
df_tsne = pd.DataFrame(sklearn_transf, columns=['tsne_x','tsne_y'])
df_tsne['cells'] = [sample[:-len('_mcc')] for sample in df.columns.tolist()]

if FLIP_AXES:
	x = df_tsne.tsne_x.values
	y = df_tsne.tsne_y.values
	cells = df_tsne.cells.values

	df_tsne = pd.DataFrame()
	df_tsne['tsne_x'] = y
	df_tsne['tsne_y'] = x
	df_tsne['cells'] = cells

df_tsne.to_csv(outfile, sep="\t", index=False)
print("Saved output to file %s" % outfile)


# plt.plot(sklearn_transf[:,0], sklearn_transf[:,1], '.r')
# plot
fig, ax = plt.subplots()
ax.scatter(df_tsne.tsne_x, df_tsne.tsne_y, s=1)
plt.show()


#!/usr/bin/env python3

"""Generate tSNE coordinates
"""
# import matplotlib
# matplotlib.use('Agg')
# import matplotlib.pyplot as plt

# import numpy as np
# import pandas as pd
# import time
# import os
# import logging
# import glob

from sklearn.decomposition import PCA
from sklearn.manifold import TSNE 

from __init__ import *
from snmcseq_utils import create_logger
from snmcseq_utils import plot_tsne_values
from snmcseq_utils import plot_tsne_labels



def run_tsne(df, perp=30, n_pc=50, n_tsne=2, 
             random_state=1, output_file=None):
    """run tsne on "_mcc$" columns
    """
    ti = time.time()
    
    df = df.filter(regex='_mcc$')
    logging.info("Running tsne: {} PC, {} perp, {} dim.\nInput shape: {}".format(n_pc, perp, n_tsne, df.shape))
    
    pca = PCA(n_components=n_pc)
    pcs = pca.fit_transform(df.T)

    tsne = TSNE(n_components=n_tsne, init='pca', random_state=random_state, perplexity=perp, verbose=0)
    ts = tsne.fit_transform(pcs)
 
    if n_tsne == 2: 
        df_tsne = pd.DataFrame(ts, columns=['tsne_x','tsne_y'])
    elif n_tsne == 3:
        df_tsne = pd.DataFrame(ts, columns=['tsne_x','tsne_y', 'tsne_z'])

    df_tsne['sample'] = [sample[:-len('_mcc')] for sample in df.columns.tolist()]
    df_tsne = df_tsne.set_index('sample')
    
    if output_file:
        df_tsne.to_csv(output_file, sep="\t", na_rep='NA', header=True, index=True)
        logging.info("Saved tsne coordinates to file. {}".format(output_file))

    tf = time.time()
    logging.info("Done with tSNE. running time: {} seconds.".format(tf - ti))
    
    return df_tsne

# def plot_tsne(df_tsne, output_file=None, s=10, **kwargs):
# 	"""plot plain tsne coordinates
# 	"""
# 	fig, ax = plt.subplots()
# 	ax.scatter(df_tsne.tsne_x, df_tsne.tsne_y, s=s, **kwargs)
# 	if output_file:
# 		title = os.path.basename(output_file).split('.')[0]
# 		ax.set_title(title)
# 		fig.savefig(output_file)
# 		logging.info("Saved tsne plot to file. {}".format(output_file))
# 		plt.close('all')
# 	return


def run_tsne_CEMBA(ens, perps=PERPLEXITIES, n_pc=N_PC, n_dim=N_DIM):
	"""
	run default tsnes for one ensemble
	"""
	ens_path = os.path.join(PATH_ENSEMBLES, ens)
	nmcc_files = sorted(glob.glob(os.path.join(ens_path, 'binc/binc_*_nmcc_{}.tsv'.format(ens)))) 

	if not os.path.isdir(os.path.join(ens_path, 'tsne')):
		os.makedirs(os.path.join(ens_path, 'tsne'))
	if not os.path.isdir(os.path.join(ens_path, 'plots')):
		os.makedirs(os.path.join(ens_path, 'plots'))

	for nmcc_file in nmcc_files:
		nmcc_basename = os.path.basename(nmcc_file) 
		df = pd.read_table(nmcc_file, dtype={'chr': object})
		for perp in perps:
			output_coords = os.path.join(ens_path, 'tsne/tsne_ndim{}_perp{}_npc{}_{}'.format(n_dim, perp, n_pc, nmcc_basename))
			output_plot = os.path.join(ens_path, 'plots/tsne_ndim{}_perp{}_npc{}_{}.pdf'.format(n_dim, perp, n_pc, nmcc_basename[:-len('.tsv')]))
			df_tsne = run_tsne(df, perp=perp, n_pc=n_pc, n_tsne=n_dim, output_file=output_coords)
			# if n_dim == 2:
			# 	plot_tsne(df_tsne, output_file=output_plot)
	return

if __name__ == '__main__':

	ti = time()
	log = create_logger()

	enss = ['Ens1', 'Ens2', 'Ens3', 'Ens4']
	for ens in enss:	
		run_tsne_CEMBA(ens)

	tf = time()
	log.info("total tSNE running time: {} second s".format(tf-ti))

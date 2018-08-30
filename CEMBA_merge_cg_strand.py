#!/usr/bin/env python3

from __init__ import *
from natsort import natsorted

import snmcseq_utils


def merge_cg_strand(df):
	"""Merge allCG table to allCG (one strand) table
	allCG table: chr, pos, strand, context, mc, c, methylation
	"""

	df['pos'] -= (df['strand'].values == '-').astype(int)

	res = df[['chr', 'pos', 'mc', 'c']].groupby(['chr', 'pos']).sum()
	res['strand'] = '+'
	res['context'] = 'CGN'
	res['methylation'] = 1
	res = res[['strand', 'context', 'mc', 'c', 'methylation']].reset_index()
	return res

if __name__ == '__main__':

	log = snmcseq_utils.create_logger()

	ens = 'Ens0'
	cluster_type = 'mCH_npc50_k30_merged'
	files = natsorted(glob.glob(os.path.join(PATH_ENSEMBLES, ens, 'allc_merged', cluster_type, 
			'allc_merged_mCG_{}_*_{}.tsv'.format(cluster_type, ens))))	
	
	# skip files with mergestrands in the name	
	for file in files:
		if 'mergestrands' in file:
			logging.info("skip file {}".format(file))
	files = [file for file in files if 'mergestrands' not in file]

			# raise ValueError("mergestrands in file already: {}".format(file))


	for i, file in enumerate(files):
		logging.info("Processing file: {} ({}/{})".format(file, i+1, len(files)))
		df = pd.read_table(file, header=None, 
						names=['chr', 'pos', 'strand', 'context', 'mc', 'c', 'methylation'], 
						dtype={'chr': object})

		res = merge_cg_strand(df)
		res.to_csv(os.path.splitext(file)[0]+'_mergestrands.tsv', sep='\t', na_rep='NA', 
					header=False, index=False)

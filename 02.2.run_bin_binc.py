#!/usr/bin/env python3

import pandas as pd
import os
import numpy as np
import multiprocessing as mp
import glob

import snmcseq_utils


def combine_binc(sample_path, output_filename=None, 
	chromosomes=None, species='human',
	input_binsize=10000, output_binsize=100000):
	"""
	combine 10kb bins to 100kb bins
	"""
	sample = os.path.basename(sample_path)

	if not output_filename:
		output_filename = (sample_path + "/binc_" + sample + "_" 
						+ str(output_binsize) + '_allchr' + ".tsv")

	if os.path.isfile(output_filename):
		print("File exists "+output_filename+", skipping...")
		return 0
	else:
		print("Processing: " + output_filename)

	if chromosomes == None:
	    if species == 'human':
	        chromosomes = snmcseq_utils.get_human_chromosomes()
	    elif species == 'mouse':
	        chromosomes = snmcseq_utils.get_mouse_chromosomes()

	mC_summary = pd.DataFrame(columns=[['chr', 'bin', 'mCG', 'CG', 'mCH', 'CH']])	
	# chrosome ordered '1,2,3,...,9,10,...,22,X'
	for chromosome in chromosomes:
		fname = sample_path + "/binc_" + sample + "_" + str(input_binsize) +  "_" + chromosome + ".tsv"
		if not os.path.isfile(fname):
			print("bin_allc: " + fname + " does not exist.")
			return False
            
		df = pd.read_table(fname)

		if species == 'human':
			bins = np.arange(0, snmcseq_utils.get_chrom_lengths_human()[chromosome], output_binsize)
		elif species == 'mouse':
			bins = np.arange(0, snmcseq_utils.get_chrom_lengths_mouse()[chromosome], output_binsize)

        # mC
		mC = df.groupby(pd.cut(df.bin, bins)).sum()[['mCG', 'CG', 'mCH', 'CH']].fillna(0)
		mC['chr'] = chromosome
		mC['bin'] = bins[:-1]

		mC_summary = pd.concat([mC_summary, mC], ignore_index=True)

	mC_summary = mC_summary[['chr', 'bin', 'mCG', 'CG', 'mCH', 'CH']]
	print(mC_summary.head())
	mC_summary.to_csv(output_filename, 
			na_rep='NA', sep="\t", header=True, index=False)

	return True


def run_combine_binc(binc_dir, species='human', bin_size=100000, nprocs=16):
	"""
	run individual combine_binc from 10kb bins to (100kb bins)
	"""

	local_dirs = glob.glob(os.path.join(binc_dir, '*_indexed')) 

	nprocs = min(nprocs, len(local_dirs))
	pool = mp.Pool(processes=nprocs)
	pool_results = [pool.apply_async(combine_binc, 
									args=(local_dir, ), 
									kwds={'species': species,
										'output_binsize': bin_size,
										}) 
					for local_dir in local_dirs]
	pool.close()
	pool.join()

	return pool_results


if __name__ == '__main__':

	binc_dir = './data/binc/binc_MB_EB'

	logger = snmcseq_utils.create_logger()	
	logger.info('Begin...')

	run_combine_binc(binc_dir, species='human', bin_size=100000, nprocs=8)

	logger.info('Done!')

	# sample_path = './data/binc/binc_MB_v1/160729_MB_v1_hs_25yr_MFG_pool_9_AD010_indexed'
	# combine_binc(sample_path)
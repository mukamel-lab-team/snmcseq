#!/usr/bin/env python3

import argparse
import os
import glob
import multiprocessing as mp

from bin_allc_files import bin_allc, get_sample_from_path

def wrap_bin_allc(local_dir, 
				out_dir=None,
				species='human',
				bin_size=10000,
				compressed=True):
	"""
	given the full path of the directory containing stores allc tables,
	output binned methylation levels.

	"""
	if not out_dir:
		try:
			out_dir = os.path.join(globals()['OUT_DIR'], os.path.basename(local_dir))
		except:
			raise ValueError('Output directory incorrect')	

	return bin_allc(get_sample_from_path(local_dir), path=local_dir,
					outpath=out_dir, species=species, bin_size=bin_size, 
					compressed=compressed)


def run_bin_allc(ALLC_DIR, nprocs=16):
	"""
	run bin_allc in parallel
	"""
	local_dirs = glob.glob(os.path.join(ALLC_DIR, '*_bismark')) 

	nprocs = min(nprocs, len(local_dirs))
	pool = mp.Pool(processes=nprocs)
	pool_results = pool.map(wrap_bin_allc, local_dirs)
	pool.close()
	pool.join()

	return pool_results

# def create_parser():
# 	"""
# 	generate parser
# 	"""
# 	parser = argparse.ArgumentParser()
# 	parser.add_argument('--binallc', help="bin allc files", action='store_true')

# 	# parser.add_argument('--all', help="going through all pipelines", action='store_true')
# 	return parser

if __name__ == '__main__':

	# parser = create_parser()
	# args = parser.parse_args()

	# if args.binallc:
	# 	pipe_bin_allc()

	ALLC_DIR = '/cndd/Public_Datasets/single_cell_methylome/allc_singlecells/hs_20161229'
	OUT_DIR = '/cndd/Public_Datasets/single_cell_methylome/binc/human_v1' 

	run_bin_allc(ALLC_DIR, nprocs=16)

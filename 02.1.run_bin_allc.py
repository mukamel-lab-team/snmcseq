#!/usr/bin/env python3

import argparse
import os
import glob
import multiprocessing as mp

from bin_allc_files import bin_allc
from snmcseq_utils import create_logger


def run_bin_allc(allc_dir, outdir, nprocs=16):
	"""
	run bin_allc in parallel
	"""
	# local_dirs = glob.glob(os.path.join(ALLC_DIR, '*_bismark')) 
	local_dirs = glob.glob(os.path.join(allc_dir, '*_indexed')) 

	if not os.path.exists(outdir):
	    os.makedirs(outdir)

	nprocs = min(nprocs, len(local_dirs))
	pool = mp.Pool(processes=nprocs)
	pool_results = [pool.apply_async(bin_allc, 
									args=(local_dir, ), 
									kwds={'outpath': os.path.join(outdir, os.path.basename(local_dir)), 
										'species': 'human',
										'bin_size': 10000,
										'compressed': True,
										}) 
					for local_dir in local_dirs]
					
	pool.close()
	pool.join()

	return pool_results


if __name__ == '__main__':

	# parser = create_parser()
	# args = parser.parse_args()


	ALLC_DIR = '/cndd/fangming/snmcseq_dev/data/allc/MB_EB'
	OUT_DIR = '/cndd/fangming/snmcseq_dev/data/binc/binc_MB_EB' 

	logger = create_logger()
	logger.info('Begin ...')

	run_bin_allc(ALLC_DIR, OUT_DIR, nprocs=16)

	logger.info('Done!')

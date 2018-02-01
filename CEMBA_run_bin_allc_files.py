#!/usr/bin/env python3

import argparse
import os
import glob
import multiprocessing as mp


from __init__ import *
from CEMBA_bin_allc_files import bin_allc
from snmcseq_utils import create_logger

# def bin_allc(allc_file, 
#     convention='CEMBA',
#     bin_size=10000, 
#     contexts=CONTEXTS,
#     chromosomes=None, 
#     species='mouse',
#     compression='gzip',
#     overwrite=False
#     ):

def run_bin_allc(allc_dirs, 
	contexts=CONTEXTS, 
	bin_size=BIN_SIZE,
	chromosomes=None,
	species='mouse',
	compression='gzip',
	convention='CEMBA',
	overwrite=False,
	nprocs=1):
	"""
	run bin_allc in parallel
	each allc_dir is from a CEMBA/Datasets/dataset_name/allc/
	"""
	logger = create_logger()

	for allc_dir in allc_dirs:
		assert os.path.isdir(allc_dir)

	allc_files = []
	for allc_dir in allc_dirs:
		allc_files += glob.glob(os.path.join(allc_dir, 'allc_*.tsv.bgz')) 

	nprocs = min(nprocs, len(allc_files))

	logger.info("""Begin run_bin_allc:{}\n
				Number of processes:{}\n
				Number of allc_files:{}\n
				""".format(allc_dirs, nprocs, len(allc_files)))
	
	pool = mp.Pool(processes=nprocs)
	pool_results = [pool.apply_async(bin_allc, 
									args=(allc_file, ), 
									kwds={'bin_size': bin_size,
										'contexts': contexts,
										'chromosomes': chromosomes,
										'species': species,
										'compression': compression,	
										'convention': convention,
										'overwrite': overwrite,
										}) 
					for allc_file in allc_files]
					
	pool.close()
	pool.join()

	logger.info('Done!')
	return pool_results

def create_parser():
    """

    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_allc_dirs", 
    	required=True,
    	nargs='+', 
    	help="list of allc file folders")
    parser.add_argument("-c", "--contexts", 
    	nargs='+', 
    	default=CONTEXTS, 
    	help="list of contexts: CH/CG/...")
    parser.add_argument("-b", "--bin_size", 
    	type=int,
		default=BIN_SIZE,
    	help="a bin size")
    parser.add_argument("-chr", "--chromosomes", 
    	nargs='+', 
    	default=None, 
    	help="list of chromosomes")
    parser.add_argument("-sp", "--species", 
		default='mouse',
    	help="species: mouse or human")
    parser.add_argument("-cp", "--compression", 
		default='gzip',
    	help="compression format of allc table (bgzip -> gzip)")
    parser.add_argument("-f", "--overwrite", 
    	action='store_true',
    	help="overwrites a file if it exists")
    parser.add_argument("-n", "--nprocs", 
    	type=int, 
    	default=1,
    	help="number of processes")

    return parser


if __name__ == '__main__':

	parser = create_parser()
	args = parser.parse_args()

	run_bin_allc(args.input_allc_dirs, 
		contexts=args.contexts, 
		bin_size=args.bin_size,
		chromosomes=args.chromosomes,
		species=args.species,
		compression=args.compression,
		convention='CEMBA',
		overwrite=args.overwrite,
		nprocs=args.nprocs)

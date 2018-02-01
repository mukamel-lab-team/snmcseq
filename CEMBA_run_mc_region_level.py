#!/usr/bin/env python3
"""
"""

import argparse
import os
import glob
import multiprocessing as mp

# from __init__ import *
from CEMBA_mc_region_level import mc_region_level_worker 
from snmcseq_utils import create_logger


# def mc_region_level_worker(allc_file, output_file, bed_file,
#     contexts=['CH', 'CG']):
	

def run_mc_region_level(allc_files, output_files, 
	bed_file, 
	contexts=['CH', 'CG'], 
	# overwrite=False,
	nprocs=1):
	"""
	run mc_gene_level in parallel
	"""
	logger = create_logger()

	# for allc_dir in allc_dirs:
	# 	assert os.path.isdir(allc_dir)

	# allc_files = []
	# for allc_dir in allc_dirs:
	# 	allc_files += glob.glob(os.path.join(allc_dir, 'allc_*.tsv.bgz')) 

	nprocs = min(nprocs, len(allc_files))

	logger.info("""Begin run_mc_region_level.\n
				Number of processes:{}\n
				Number of allc_files:{}\n
				""".format(nprocs, len(allc_files)))
	
	pool = mp.Pool(processes=nprocs)
	pool_results = [pool.apply_async(mc_region_level_worker, 
									args=(allc_file, output_file, bed_file), 
									kwds={'contexts': contexts, 
										# 'genebody': genebody,
										# 'convention': convention,
										# 'overwrite': overwrite,
										}) 
					for allc_file, output_file in zip(allc_files, output_files)]
					
	pool.close()
	pool.join()

	logger.info('Done!')

	return pool_results


def create_parser():
    """

    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_allc_files", 
    	required=True,
    	nargs='+', 
    	help="list of allc files")
    parser.add_argument("-o", "--output_files", 
    	required=True,
    	nargs='+', 
    	help="list of output_files")
    parser.add_argument("-b", "--bed_file", 
    	required=True,
    	help="bed file")
    parser.add_argument("-c", "--contexts", 
    	nargs='+', 
    	default=['CH', 'CG'], 
    	help="list of contexts: CH/CG/...")
    # parser.add_argument("-f", "--overwrite", 
    # 	action='store_true',
    # 	help="overwrites a file if it exists")
    parser.add_argument("-n", "--nprocs", 
    	type=int, 
    	default=1,
    	help="number of processes")

    return parser



if __name__ == '__main__':

	parser = create_parser()
	args = parser.parse_args()

	run_mc_region_level(args.input_allc_files, args.output_files, args.bed_file,
		contexts=args.contexts,
		# overwrite=args.overwrite,
		nprocs=args.nprocs)


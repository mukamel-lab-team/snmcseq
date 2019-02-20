#!/usr/bin/env python3
"""
"""

import argparse
import os
import glob
import multiprocessing as mp

from __init__ import *
from CEMBA_mc_gene_level import mc_gene_level 
from snmcseq_utils import create_logger


# def mc_gene_level(allc_file,
#     contexts=['CH', 'CG', 'CA'], 
#     genebody='/cndd/Public_Datasets/CEMBA/References/Annotation/gencode.vM16.annotation_genes.tsv',
#     convention='CEMBA', 
#     overwrite=False):

def run_mc_gene_level(allc_dirs, 
	contexts=CONTEXTS, 
	genebody=GENEBODY, 
	convention='CEMBA',
	overwrite=False,
	nprocs=1):
	"""
	run mc_gene_level in parallel
	"""
	logger = create_logger()

	for allc_dir in allc_dirs:
		if not os.path.isdir(allc_dir):
			raise ValueError("{} is not a valid directory!".format(allc_dir))

	allc_files = []
	for allc_dir in allc_dirs:
		allc_files += glob.glob(os.path.join(allc_dir, 'allc_*.tsv.bgz')) 

	nprocs = min(nprocs, len(allc_files))

	logger.info("""Begin run_mc_gene_level:{}\n
				Number of processes:{}\n
				Number of allc_files:{}\n
				""".format(allc_dirs, nprocs, len(allc_files)))
	
	pool = mp.Pool(processes=nprocs)
	pool_results = [pool.apply_async(mc_gene_level, 
									args=(allc_file, ), 
									kwds={'contexts': contexts, 
										'genebody': genebody,
										# 'convention': convention,
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
    parser.add_argument("-g", "--genebody", 
		default=GENEBODY,
    	help="file of gene body")
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

	run_mc_gene_level(args.input_allc_dirs, 
		contexts=args.contexts,
		genebody=args.genebody, 
		convention='CEMBA',
		overwrite=args.overwrite,
		nprocs=args.nprocs)


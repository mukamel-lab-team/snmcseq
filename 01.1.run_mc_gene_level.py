#!/usr/bin/env python3
"""
Usage:

1. make sure allc_path/sample_name ends with "_indexed"
2. context should be in consistent with output file names
"""

import argparse
import os
import glob
import multiprocessing as mp

from mc_gene_level import mc_gene_level 
from snmcseq_utils import create_logger

# def mc_gene_level(sample,
#     genebody='/cndd/projects/Public_Datasets/references/hg19/transcriptome/gencode.v19.annotation_genes_mypy.tsv',
#     outdir="./genebody",
#     context='CH')

def run_mc_gene_level(allc_dir, outdir, context, 
	genebody='/cndd/projects/Public_Datasets/references/hg19/transcriptome/gencode.v19.annotation_genes_mypy.tsv', 
	nprocs=16):
	"""
	run mc_gene_level in parallel
	"""
	# local_dirs = glob.glob(os.path.join(allc_dir, '*_bismark')) 
	local_dirs = glob.glob(os.path.join(allc_dir, '*_indexed')) 

	if not os.path.exists(outdir):
	    os.makedirs(outdir)

	nprocs = min(nprocs, len(local_dirs))
	pool = mp.Pool(processes=nprocs)
	pool_results = [pool.apply_async(mc_gene_level, 
									args=(local_dir, ), 
									kwds={'outdir': outdir, 
										'context': context,
										'genebody': genebody}) 
					for local_dir in local_dirs]
					
	pool.close()
	pool.join()

	return pool_results


if __name__ == '__main__':

	# parser = create_parser()
	# args = parser.parse_args()

	CONTEXT = 'CH'
	ALLC_DIR = './data/allc/MB_EA'
	OUT_DIR = './data/gene_level/genebody_m%s_MB_EA' % CONTEXT 


	logger = create_logger()
	logger.info('Begin ...')

	run_mc_gene_level(ALLC_DIR, OUT_DIR, CONTEXT, nprocs=16)

	logger.info('Done!')

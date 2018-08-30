#!/usr/bin/env python3

import argparse
import os
import glob
import multiprocessing as mp

from mc_gene_level import mc_gene_level 


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
	local_dirs = glob.glob(os.path.join(allc_dir, '*_bismark')) 

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


	ALLC_DIR = '/cndd/Public_Datasets/single_cell_methylome/allc_singlecells/hs_MB_EB'
	OUT_DIR = '/cndd/Public_Datasets/single_cell_methylome/gene_level/human/genebody_CG_hv2' 
	CONTEXT = 'CG'

	run_mc_gene_level(ALLC_DIR, OUT_DIR, CONTEXT, nprocs=16)

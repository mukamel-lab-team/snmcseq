#!/usr/bin/env python3
"""
"""

from __init__ import *
import multiprocessing as mp
import argparse

import snmcseq_utils
from CEMBA_mc_region_level import mc_region_level_worker 
from snmcseq_utils import create_logger


# def mc_region_level_worker(allc_file, output_file, bed_file,
#     contexts=['CH', 'CG']):
	

def run_mc_region_level(allc_files, output_files, 
	bed_file, 
	contexts=CONTEXTS,
	compress=True, 
	cap=2,
	# overwrite=False,
	nprocs=1):
	"""
	run mc_gene_level in parallel

	cap: remove totalC>cap
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
				Bed file: {}\n
				""".format(nprocs, len(allc_files), bed_file))
	
	pool = mp.Pool(processes=nprocs)
	pool_results = [pool.apply_async(mc_region_level_worker, 
									args=(allc_file, output_file, bed_file), 
									kwds={'contexts': contexts, 
										'compress': compress,
										'cap': cap,
										}) 
					for allc_file, output_file in zip(allc_files, output_files)]
					

	pool.close()
	pool.join()
	# mc_region_level_worker(allc_files[0], output_files[0], bed_file, contexts=contexts, compress=compress)	


	return pool_results


def run_mc_region_level_CEMBA(dataset, bed_file, output_dirname, contexts=CONTEXTS, compress=True, nprocs=1, cap=2):
	"""Generate bed_file and 
	"""
	assert snmcseq_utils.isdataset(dataset)

	if not os.path.isdir(os.path.join(PATH_DATASETS, dataset, output_dirname)):
		logging.info("Creating directory: {}".format(os.path.join(PATH_DATASETS, dataset, output_dirname)))
		os.makedirs(os.path.join(PATH_DATASETS, dataset, output_dirname))


	allc_files = sorted(glob.glob(os.path.join(PATH_DATASETS, dataset, 'allc', 'allc_*.tsv.bgz')))
	cells = [os.path.basename(allc_file)[len('allc_'):-len('.tsv.bgz')] for allc_file in allc_files]
	output_files = [os.path.join(PATH_DATASETS, dataset, output_dirname, cell+'.tsv') for cell in cells]

	logging.info("mc region level counts: {}".format(dataset))



	run_mc_region_level(allc_files, output_files, 
		bed_file, 
		contexts=CONTEXTS,
		compress=compress, 
		nprocs=nprocs, 
		cap=cap,
	)

	return



def create_parser():
    """

    """
    parser = argparse.ArgumentParser()
    # parser.add_argument("-i", "--input_allc_files", 
    # 	required=True,
    # 	nargs='+', 
    # 	help="list of allc files")
    # parser.add_argument("-o", "--output_files", 
    # 	required=True,
    # 	nargs='+', 
    # 	help="list of output_files")
    parser.add_argument("-d", "--datasets", 
    	required=True,
    	nargs='+', 
    	help="list of datasets")
    parser.add_argument("-o", "--output_dirname", 
    	required=True,
    	help="output dirname UNDER dataset folders (example: gene_level, binc, atac_peaks, ...)")
    parser.add_argument("-b", "--bed_file", 
    	required=True,
    	help="bed file")
    parser.add_argument("-c", "--contexts", 
    	nargs='+', 
    	default=CONTEXTS, 
    	help="list of contexts: CH/CG/... default: CH CG CA")

    parser.add_argument("-cp", "--cap", 
    	type=int, 
    	default=2,
    	help="Exclude totalC greater than this values (default 2)")

    parser.add_argument("-n", "--nprocs", 
    	type=int, 
    	default=1,
    	help="number of processes (default 1)")

    return parser



if __name__ == '__main__':

	log = create_logger()
	parser = create_parser()
	args = parser.parse_args()

	datasets = args.datasets
	bed_file = args.bed_file
	output_dirname = args.output_dirname
	contexts = args.contexts
	nprocs = args.nprocs
	cap = args.cap


	if '/' in output_dirname:
		raise ValueError("'/' is not allowed in a directory name, choose names such as: atac_peak, dmr, etc...")
	elif 'binc' in output_dirname:
		raise ValueError("'binc' is not allowed as a directory name!")
	elif 'gene_level' in output_dirname:
		raise ValueError("'gene_level' is not allowed as a directory name!")

	# datasets = ['CEMBA_3C_171206', 'CEMBA_3C_171207', 'CEMBA_4B_171212', 'CEMBA_4B_171213', 'CEMBA_4B_180104']	
	# bed_files = ['/cndd/Public_Datasets/CEMBA/snATACSeq/Datasets/{}/{}_merged.bed'.format(dataset, dataset) for dataset in datasets]
	# output_dirname = 'atac_peak_regions'
	# nprocs = 16

	logging.info(""" CEMBA mc region level counting:
			Datasets: {}
			Bed file: {}
			Output dirname: {}
			Cap: {}
			Number of processes: {}
		""".format(datasets, bed_file, output_dirname, cap, nprocs))

	for dataset in datasets:
		assert snmcseq_utils.isdataset(dataset)
	assert os.path.isfile(bed_file)

	for dataset in datasets:
		res = run_mc_region_level_CEMBA(dataset, bed_file, output_dirname, 
			contexts=CONTEXTS, compress=True, cap=cap, nprocs=nprocs)

	# run_mc_region_level(args.input_allc_files, args.output_files, args.bed_file,
	# 	contexts=args.contexts,
	# 	# overwrite=args.overwrite,
	# 	nprocs=args.nprocs)

#!/usr/bin/env python3
"""Organize files before running the pipeline 
"""

from __init__ import *
import multiprocessing as mp

from snmcseq_utils import create_logger
from snmcseq_utils import cd
from snmcseq_utils import isdataset 


def gzip_to_bgzip(src, dst):
	"""src is a gzip file, make it bgzip at dst
	"""
	logging.info("Begin gzip to bgzip: {} -> {}".format(src, dst))

	# .tsv.gz to .tsv 
	cmd = 'gzip -cd {} > {}'.format(src, os.path.splitext(src)[0]) 
	os.system(cmd)

	# bgzip
	cmd = 'bgzip -c {} > {}'.format(os.path.splitext(src)[0], dst) 
	os.system(cmd)

	# rm tmp files
	os.remove(os.path.splitext(src)[0])
	logging.info("Done gzip to bgzip: {} -> {}".format(src, dst))

	return

def CEMBA_gzip_to_bgzip(dataset, nprocs=1):
	"""Given a dataset, get paired allc src (under from_ecker_lab folder) 
	and dst (under allc folder) list, and send it to gzip_to_bgzip
	"""
	dataset_folder = os.path.join(PATH_DATASETS, dataset, 'from_ecker_lab')	
	assert os.path.isdir(dataset_folder)

	with cd(dataset_folder):

		# get srcs
		possible_allc_files = glob.glob('allc_*.tsv.gz')
		if possible_allc_files:
			srcs = possible_allc_files
		else:
			allc_folders = [item for item in glob.glob('*') if os.path.isdir(item)] 
			srcs = []
			for allc_folder in allc_folders:
				possible_allc_files = glob.glob('./{}/allc_*.tsv.gz'.format(allc_folder))
				assert len(possible_allc_files) == 1
				srcs.append(possible_allc_files[0])

		# get dsts based on srcs
		if not os.path.isdir(os.path.join(PATH_DATASETS, dataset, 'allc')):
			os.makedirs(os.path.join(PATH_DATASETS, dataset, 'allc'))
		dsts = [os.path.join(PATH_DATASETS, dataset, 'allc', os.path.splitext(os.path.basename(src))[0]+'.bgz')
				for src in srcs]
		n_files = len(srcs)

		# do gzip to bgzip
		# for i, (src, dst) in enumerate(zip(srcs, dsts)): # can be parallelized
		# 	logging.info("Begin Gzip to bgzip: {} -> {} ({}/{})".format(src, dst, i+1, n_files))
		# 	gzip_to_bgzip(src, dst)
		# 	logging.info("Done Gzip to bgzip: {} -> {} ({}/{})".format(src, dst, i+1, n_files))

		# do gzip to bgzip (parallelized)
		nprocs = min(nprocs, len(srcs))
		logging.info("""Begin organize gzip to bgzip\n
					Number of processes:{}\n
					Number of allc_files:{}\n
					""".format(nprocs, len(srcs)))
		
		pool = mp.Pool(processes=nprocs)
		pool_results = [pool.apply_async(gzip_to_bgzip, 
										args=(src, dst), 
										) 
						for (src, dst) in zip(srcs, dsts)]
						
		pool.close()
		pool.join()

		logging.info('Done!')
	return

def tabix_index_allc(src, skip_lines=0):
	"""src is a .bgz allc table
	"""
	os.system('tabix -f -s 1 -b 2 -e 2 -S {} {}'.format(skip_lines, src)) 
	return 

def main(dataset, nprocs=1):
	# make allc folder and get files from_ecker_lab
	if not os.path.isdir(os.path.join(PATH_DATASETS, dataset, 'allc')):
		os.makedirs(os.path.join(PATH_DATASETS, dataset, 'allc'))

	# gzip to bgzip and have files in place
	CEMBA_gzip_to_bgzip(dataset, nprocs=nprocs)

	# tabix index
	with cd(os.path.join(PATH_DATASETS, dataset, 'allc')):
		allc_files = glob.glob(os.path.join(PATH_DATASETS, dataset, 'allc', 'allc_*.tsv.bgz'))
		for i, allc_file in enumerate(allc_files):	# try tabix (could be parallelized)
			logging.info("Begin tabix indexing ({}/{}): {}".format(i+1, len(allc_files), allc_file))
			tabix_index_allc(allc_file)
			logging.info("Done tabix indexing ({}/{}): {}".format(i+1, len(allc_files), allc_file))

if __name__ == '__main__':

	log = create_logger()
	
	datasets = [
		'CEMBA_2A_180122',
		'CEMBA_2A_180123',
		'CEMBA_2C_180409',
		'CEMBA_2C_180410',
		'CEMBA_3B_180312',
		'CEMBA_3B_180501',
		'CEMBA_3D_180412',
		'CEMBA_3D_180416',
		'CEMBA_5B_180514',
		'CEMBA_5B_180529',
		'CEMBA_5D_180605',
		'CEMBA_5D_180612',
		'CEMBA_7B_180423',
		'CEMBA_7B_180424',
	]

	nprocs = 8 
	# check if dataset exists 
	for dataset in datasets:
		if not isdataset(dataset):
			raise ValueError('Dataset {} not found in datasets'.format(dataset))

	# process the dataset (organize files, gzip to bgzip, and tabix index)
	for dataset in datasets:
		logging.info("Begin organizing files for dataset: {}".format(dataset))
		main(dataset, nprocs)
		logging.info("Done organizing files for dataset: {}".format(dataset))



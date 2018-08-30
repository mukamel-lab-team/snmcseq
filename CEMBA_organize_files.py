#!/usr/bin/env python3
"""Organize files before running the pipeline 
"""

from __init__ import *
from snmcseq_utils import create_logger
from snmcseq_utils import cd
from snmcseq_utils import isdataset 


def gzip_to_bgzip(src, dst):
	"""src is a gzip file, make it bgzip at dst
	"""

	# .tsv.gz to .tsv 
	cmd = 'gzip -cd {} > {}'.format(src, os.path.splitext(src)[0]) 
	os.system(cmd)

	# bgzip
	cmd = 'bgzip -c {} > {}'.format(os.path.splitext(src)[0], dst) 
	os.system(cmd)

	# rm tmp files
	os.remove(os.path.splitext(src)[0])

	return

def CEMBA_gzip_to_bgzip(dataset):
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
		for i, (src, dst) in enumerate(zip(srcs, dsts)): # can be parallelized
			logging.info("Begin Gzip to bgzip: {} -> {} ({}/{})".format(src, dst, i+1, n_files))
			gzip_to_bgzip(src, dst)
			logging.info("Done Gzip to bgzip: {} -> {} ({}/{})".format(src, dst, i+1, n_files))
	return

def tabix_index_allc(src, skip_lines=0):
	"""src is a .bgz allc table
	"""
	os.system('tabix -f -s 1 -b 2 -e 2 -S {} {}'.format(skip_lines, src)) 
	return 

def main(dataset):
	# make allc folder and get files from_ecker_lab
	if not os.path.isdir(os.path.join(PATH_DATASETS, dataset, 'allc')):
		os.makedirs(os.path.join(PATH_DATASETS, dataset, 'allc'))

	# gzip to bgzip and have files in place
	CEMBA_gzip_to_bgzip(dataset)

	# tabix index
	with cd(os.path.join(PATH_DATASETS, dataset, 'allc')):
		allc_files = glob.glob(os.path.join(PATH_DATASETS, dataset, 'allc', 'allc_*.tsv.bgz'))
		for i, allc_file in enumerate(allc_files):	# try tabix (could be parallelized)
			logging.info("Begin tabix indexing ({}/{}): {}".format(i+1, len(allc_files), allc_file))
			tabix_index_allc(allc_file)
			logging.info("Done tabix indexing ({}/{}): {}".format(i+1, len(allc_files), allc_file))

if __name__ == '__main__':

	log = create_logger()
	datasets = ['CEMBA_3F_180109']

	# check if dataset exists 
	for dataset in datasets:
		if not isdataset(dataset):
			raise ValueError('Dataset {} not found in datasets'.format(dataset))

	# process the dataset (organize files, gzip to bgzip, and tabix index)
	for dataset in datasets:
		logging.info("Begin organizing files for dataset: {}".format(dataset))
		main(dataset)
		logging.info("Done organizing files for dataset: {}".format(dataset))



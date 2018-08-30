#!/usr/bin/env python3

import subprocess as sp

from __init__ import *
from snmcseq_utils import create_logger


def call_DMR_wrapper(allc_paths, output_prefix, nprocs=8):
	"""
	"""	
	samples = [os.path.basename(allc_path)[len('allc_'):-len('.tsv')] for allc_path in allc_paths]
	# nothing should be there (including space) after each line in cmd
	cmd = (
	"""methylpy DMRfind \\
		--allc-files {} \\
		--samples {} \\
		--output-prefix {} \\
		--mc-type "CGN" \\
		--num-procs {} \\
		--min-num-dms 0 \\
		--dmr-max-dist 250""".format(' '.join(allc_paths), ' '.join(samples), output_prefix, nprocs)
	)

	logging.info("*"*10)	
	logging.info("\n"+cmd)	
	logging.info("*"*10)	
	try:
		sp.run(cmd, shell=True)
	except:
		sp.call(cmd, shell=True)
	return

if __name__ == '__main__':

	log = create_logger()


	cluster_type = 'mCH_npc50_k30_merged'
	ens = 'Ens0'
	merge_strands = True

	if merge_strands:
		merge_strands = '_mergestrands'
	else:
		merge_strands = ''


	ens_path = os.path.join(PATH_ENSEMBLES, ens)
	allc_paths = sorted(glob.glob(os.path.join(ens_path, 'allc_merged', 
						cluster_type, 'allc_merged_mCG_{}_*_{}{}.tsv'.format(cluster_type, ens, merge_strands))))

	output_path = os.path.join(ens_path, 'dmr')
	output_prefix = os.path.join(output_path, 'dmr_allc_merged_mCG_{}'.format(cluster_type))
	if not os.path.isdir(output_path):
		os.makedirs(output_path)
		logging.info("Created path: {}".format(output_path))

	call_DMR_wrapper(allc_paths, output_prefix, nprocs=8)









	# cluster_type = 'cluster_mCHmCG_lv_npc50_k5'
	# enss = ['Ens10']
	# for ens in enss: 
	# 	ens_path = os.path.join(PATH_ENSEMBLES, ens)
	# 	allc_paths = sorted(glob.glob(os.path.join(ens_path, 'allc_merged/allc_merged_mCG_{}_*_{}.tsv'.format(cluster_type, ens))))

	# 	output_path = os.path.join(ens_path, 'dmr')
	# 	output_prefix = os.path.join(output_path, 'dmr_allc_merged_mCG_{}'.format(cluster_type))
	# 	if not os.path.isdir(output_path):
	# 		os.makedirs(output_path)
	# 		logging.info("Created path: {}".format(output_path))

	# 	call_DMR_wrapper(allc_paths, output_prefix, nprocs=8)



#!/usr/bin/env python3

from __init__ import *

from natsort import natsorted

from snmcseq_utils import create_logger
import CEMBA_merge_allc
import CEMBA_call_DMR
import CEMBA_allc_dmr_to_mc_single
import CEMBA_setup_annoj_v2

if __name__ == '__main__':

	log = create_logger()

	enss = ['Ens35']  # example
	nprocs = 4
	# choose a clustering type
	cluster_type = 'cluster_mCHmCG_lv_npc50_k30'
	MERGE_ALLC = True
	CALL_DMR = True
	UPLOAD_TO_MYSQL = True
	GEN_ANNOJ_TRACKS = True

	for ens in enss: 
		# ens = 'Ens1'
		# merge allc
		if MERGE_ALLC:
			CEMBA_merge_allc.merge_allc_CEMBA(ens, context='CG', cluster_type=cluster_type, database=DATABASE, 
				nprocs=nprocs,
				chunksize=1000000, n_chunk1=50, n_chunk2=20)

		# call DMRs
		if CALL_DMR:
			allc_paths = natsorted(glob.glob(os.path.join(PATH_ENSEMBLES, ens, 
							'allc_merged/allc_merged_mCG_{}_*_{}.tsv.gz'.format(cluster_type, ens))))
			output_path = os.path.join(PATH_ENSEMBLES, ens, 'dmr')
			output_prefix = os.path.join(output_path, 'dmr_allc_merged_mCG_{}'.format(cluster_type))

			if not os.path.isdir(output_path):
				os.makedirs(output_path)
				logging.info("Created path: {}".format(output_path))
			CEMBA_call_DMR.call_DMR_wrapper(allc_paths, output_prefix, nprocs=nprocs)

			# generate bed files for each cluster

		# upload data to mysql
		if UPLOAD_TO_MYSQL:
			# mc dmr to mc_single format
			allc_files = allc_paths			
			use_dmrs = True
		    dmr_files = natsorted(glob.glob(os.path.join(
		    	PATH_ENSEMBLES, 'dmr', 'cgdmr_multimodal_v2/multimodal_v2_*.bed')))  #
		    output_dir = os.path.join(PATH_ENSEMBLES, ens, 'annoj') 
		    if not os.path.isdir(output_dir):
		    	os.path.makedirs(output_dir)
		    allc_dmr_to_mc_single(allc_files, dmr_files, output_dir, use_dmrs)

		    # upload to mysql
		    # 02.mc_single_load_mysql



		# gen browser tracks
		if GEN_ANNOJ_TRACKS:
			# setup annoj


			pass

		# transfer data to brainome

		break

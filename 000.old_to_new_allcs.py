#!/usr/bin/env python3
"""Old allc format to new allc format
"""

from __init__ import *
from natsort import natsorted

import snmcseq_utils
import subprocess
import CEMBA_organize_files


if __name__ == '__main__':

	log = snmcseq_utils.create_logger()

	datasets = ['MB_v1', 'MB_EA', 'MB_EB']
	for dataset in datasets:
		logging.info("Transfer dataset: {}".format(dataset))
		orig_allcdir = '/cndd/fangming/snmcseq_dev/data/allc/{}'.format(dataset)
		path_dst = '/cndd/Public_Datasets/human_snmcseq/Datasets/{}'.format(dataset)

		with snmcseq_utils.cd(orig_allcdir):
			allc_dirs = sorted(glob.glob("*")) 
			for i, allc_dir in enumerate(allc_dirs):
				logging.info("Progress: {}/{}".format(i+1, len(allc_dirs)))
				# allc_dir is cell name
				dst = os.path.join(path_dst, 'allc', 'allc_'+allc_dir+'.tsv')
				# do not overwrite
				if os.path.isfile(dst) or os.path.isfile(dst+'.bgz'):
					continue

				with snmcseq_utils.cd(allc_dir):
					allc_files = natsorted(glob.glob("allc_*.tsv.gz"))
					assert len(allc_files) == 24

					for allc_file in allc_files:
						# decompress, remove first line, and  append
						cmd = 'gzip -cd {} | tail -n +2 >> {}'.format(allc_file, dst)
						os.system(cmd)

					# compress
					cmd = 'bgzip -c {} > {}'.format(dst, dst+'.bgz') 
					os.system(cmd)

					# remove
					cmd = 'rm {}'.format(dst)
					os.system(cmd)

					# tabix
					CEMBA_organize_files.tabix_index_allc(dst+'.bgz', skip_lines=0)

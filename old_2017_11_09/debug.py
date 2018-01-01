#!/usr/bin/env python3

import os
import glob

# DIR = '/cndd/Public_Datasets/single_cell_methylome/allc_singlecells/hs_20161229'
# for allc_dir in glob.glob(os.path.join(DIR, '*_bismark')):
# 	print('.', end='', flush=True)
# 	if not os.listdir(allc_dir):
# 		print("Empty dir: %s" % os.path.basename(allc_dir))
# 		raise

# print('Done!')


DIR = '/cndd/Public_Datasets/single_cell_methylome/binc/human_v1'

for binc_dir in glob.glob(os.path.join(DIR, 'Pool_*')):
	print('.', end='', flush=True)
	if not os.listdir(binc_dir):
		print("Empty dir: %s" % os.path.basename(binc_dir))
		

print('Done!')

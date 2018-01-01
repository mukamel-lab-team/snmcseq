#!/usr/bin/env python3

import os
# import shutil


### step 1. symbolic link and add prefix
# src_dirs = [
# 	'/cndd/Public_Datasets/single_cell_methylome/gene_level/human/genebody_mCG_hs1',
# 	'/cndd/Public_Datasets/single_cell_methylome/gene_level/human/genebody_mCG_hs2',
# 	]

src_dirs = [
	'/cndd/Public_Datasets/single_cell_methylome/binc/human_v1',
	'/cndd/Public_Datasets/single_cell_methylome/binc/human_MB_EB',
	]
biosamples = ['hv1', 
			 'hv2']

dst_dir = ('/cndd/Public_Datasets/single_cell_methylome/binc/human_combined')

if not os.path.exists(dst_dir):
	os.makedirs(dst_dir)

for bs, src_dir in zip(biosamples, src_dirs): # biosamples/src dirs
	for file in os.listdir(src_dir):	# files in src dir
		# create symbolic link over with a prefix
		os.symlink(os.path.join(src_dir, file),
				os.path.join(dst_dir, bs+'_'+file))

### step 2. merge metadata 
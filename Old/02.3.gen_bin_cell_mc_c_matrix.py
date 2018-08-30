#!/usr/bin/env python3

"""
extract and combine bin mC info of each cells 
to a matrix of gene*cells with the following format:

chr bin sample1_mc sample1_c sample2_mc sample2_c ....

"""

import pandas as pd
import os
import glob

from snmcseq_utils import create_logger

def gen_bin_cell_mc_c_matrix(binc_dir, output_fname=None, bin_size=100000, context='CH'):
	"""
	"""

	for i, sample_dir in enumerate(sorted(os.listdir(binc_dir))):

		sample_name = sample_dir 
		fname = ('binc_' + sample_name + '_' + str(bin_size) + '_' + 'allchr.tsv')
		f_fullpath = os.path.join(os.path.join(binc_dir, sample_dir), fname)

		# multi-index works similar to single index
		df = pd.read_table(f_fullpath, header=0, index_col=['chr', 'bin'])	
		# df.sort_index(inplace=True)
		# first file: generate a new dataframe that stores everything
		if i == 0:
			df_new = pd.DataFrame(index=df.index)
			df_new_cols = []

		assert df_new.index.tolist() == df.index.tolist()
		df_new[sample_name+'_mc'] = df['m'+context] 
		df_new[sample_name+'_c'] = df[context] 
		df_new_cols.append(sample_name+'_mc')
		df_new_cols.append(sample_name+'_c')
		print('Loaded sample %s.' % sample_name)

	df_new = df_new[df_new_cols]
	print(df_new.shape)

	if output_fname:
		df_new.to_csv(output_fname, 
				sep='\t', na_rep='NA', index=True, header=True)
		print('Saved to %s.' % output_fname)

	return df_new


if __name__ == '__main__':

	logger = create_logger()
	logger.info('Begin ...')

	context = 'CH'
	biosample = 'MB_EB'
	bin_size = 100000
	binc_dir = './data/binc/binc_%s' % (biosample)
	output_fname = './data/binc/binc_m%s_%s_%s_summary.tsv' % (context, biosample, str(bin_size))

	gen_bin_cell_mc_c_matrix(binc_dir, 
							output_fname=output_fname, 
							bin_size=bin_size,
							context=context)

	logger.info('Done!')



# #for test purpose
# df_2 = pd.read_table(output_fname, header=0, index_col=0)
# print(df_2.head())
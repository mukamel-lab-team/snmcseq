#!/usr/bin/env python3

"""
extract and combine genebody mC info of each cells 
to a matrix of gene*cells with the following format:

id sample1_mc sample1_c sample2_mc sample2_c ....

"""

import pandas as pd
import os
import glob


def gen_gene_cell_mc_c_matrix(genebody_dir, output_fname=None, context='CH'):
	"""
	"""

	for i, fname in enumerate(sorted(os.listdir(genebody_dir))):
		assert fname.endswith('_m%s_genebody.txt' % context)
		sample_name = fname[:-len('_m%s_genebody.txt' % context)]

		f_fullpath = os.path.join(genebody_dir, fname)
		df = pd.read_table(f_fullpath, header=0, index_col='id')	
		df.sort_index(inplace=True)
		# first file: generate a new dataframe that stores everything
		if i == 0:
			df_new = pd.DataFrame(index=df.index)
			df_new_cols = []

		assert df_new.index.tolist() == df.index.tolist()
		df_new[sample_name+'_mc'] = df['mc'] 
		df_new[sample_name+'_c'] = df['c'] 
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
	genebody_dir = '/cndd/Public_Datasets/single_cell_methylome/gene_level/human/genebody_combined'
	output_fname = '/cndd/Public_Datasets/single_cell_methylome/gene_level/human/genebody_mCH_combined-hs1-hs2.mat'
	context = 'CH'
	gen_gene_cell_mc_c_matrix(genebody_dir, output_fname=output_fname, context=context)
	print('Done!')


# #for test purpose
# df_2 = pd.read_table(output_fname, header=0, index_col=0)
# print(df_2.head())
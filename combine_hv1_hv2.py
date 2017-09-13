#!/usr/bin/env python3

"""
generate combined metadata and binned mCH file

use 'hv1' and 'hv2' header to distinguish two different batch
"""

import pandas as pd
import numpy as np

def combine_metadata(meta_hv1_fname, meta_hv2_fname, output_fname):
	"""
	"""
	df_hv1 = pd.read_table(meta_hv1_fname, header=0)	
	#####
	# df_hv1['Sample'] = ['hv1_'+item+'_R1' for item in df_hv1['Sample'].tolist()]
	df_hv1['Sample'] = ['hv1_'+item for item in df_hv1['Sample'].tolist()]
	df_hv1['Batch'] = ['hv1']*df_hv1.shape[0]
	# df_hv1 = df_hv1[['Sample', 'mCH/CH', 'Batch']]
	# print(df_hv1.head())
	# print(df_hv1.shape)
	
	df_hv2 = pd.read_table(meta_hv2_fname, header=0)	
	df_hv2['Sample'] = ['hv2_'+item for item in df_hv2['Sample'].tolist()]
	df_hv2['Batch'] = ['hv2']*df_hv2.shape[0]
	# df_hv2 = df_hv2[['Sample', 'mCH/CH', 'Batch']]
	# print(df_hv2.head())
	# print(df_hv2.shape)

	df_combined = pd.concat([df_hv1, df_hv2], ignore_index=True)
	# print(df_combined.head())
	# print(df_combined.shape)

	df_combined.to_csv(output_fname, sep='\t', header=True, index=False)

	return


def combine_bin_mCH(binc_hv1_fname, binc_hv2_fname, output_fname):	
	"""
	"""
	df_hv1 = pd.read_table(binc_hv1_fname, header=0, dtype={'chr': np.unicode, 'bin': np.int32})	
	# chr and bin not affected
	df_hv1.columns = ['chr', 'bin'] + ['hv1_'+item for item in df_hv1.columns.tolist()[2:]] 

	df_hv2 = pd.read_table(binc_hv2_fname, header=0, dtype={'chr': np.unicode, 'bin': np.int32})	
	# chr and bin not affected
	df_hv2.columns = ['chr', 'bin'] + ['hv2_'+item for item in df_hv2.columns.tolist()[2:]] 

	assert df_hv1['chr'].tolist() == df_hv2['chr'].tolist()
	assert df_hv1['bin'].tolist() == df_hv2['bin'].tolist()

	df_combined = pd.concat([df_hv1, df_hv2.iloc[:, 2:]], axis=1) 
	print(df_hv1.shape)
	print(df_hv2.shape)
	print(df_combined.shape)
	df_combined.to_csv(output_fname, sep='\t', header=True, index=False)

	return

if __name__ == '__main__':
	meta_hv1_fname = './metadata/human_v1_metadata_cells.tsv'
	meta_hv2_fname = './metadata/MB_EB_metadata_cells.tsv'
	meta_output_fname = './metadata/human_hv1_hv2_metadata_cells.tsv' 
	combine_metadata(meta_hv1_fname, meta_hv2_fname, meta_output_fname)

	# binc_hv1_fname = './tsne/human_combined_v1_CH_100000.tsv'
	# binc_hv2_fname = './tsne/human_combined_MB_EB_CH_100000.tsv'
	# binc_output_fname = './tsne/human_combined_hv1_hv2_CH_100000.tsv' 
	# combine_bin_mCH(binc_hv1_fname, binc_hv2_fname, binc_output_fname)

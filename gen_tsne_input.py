#!/usr/bin/env python3
import os
import glob
import pandas as pd

DIR = '/cndd/fangming/side_project/snmcseq_dev/temp_data'
fnames = glob.glob(os.path.join(DIR, 'binc_Pool_2256_*_1.tsv'))


columns = []
df_combined = pd.DataFrame()

for fname in fnames:
	sample_name = os.path.basename(fname).split('.')[0]
	df = pd.read_table(fname, header=None)
	df = df[[4,5]].head(100)

	df_combined[sample_name+'_mc'] = df[4].values	
	df_combined[sample_name+'_c'] = df[5].values 

# output format:
# columns: sample1_mc, sample1_c, sample2_mc, sample2_c, ...
df_combined.to_csv('test_tsne.tsv', sep='\t', index=False)

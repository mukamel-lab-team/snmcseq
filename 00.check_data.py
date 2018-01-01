#!/usr/bin/env python3

"""
1. check if the metadata and allc files are in place
2. compare their sample list
3. give stats about samples 
"""

import subprocess as sp
import pandas as pd
import os

# # rename files
# df_re = pd.read_table('/cndd/chongyuan/MB_samples_rename.txt', header=None)
# print(df_re.shape)

# path1 = './data/allc/MB_v1'
# path2 = './data/allc/MB_EA'
# path3 = './data/allc/MB_EB'

# samples1 = os.listdir(path1)
# samples2 = os.listdir(path2)
# samples3 = os.listdir(path3)







# # match metadata
# df_m1 = pd.read_table('./data/metadata/metadata_MB_v1_updated.tsv')
# df_m2 = pd.read_table('./data/metadata/metadata_MB_EA_updated.tsv')
# df_m3 = pd.read_table('./data/metadata/metadata_MB_EB_updated.tsv')

# print(samples3[:3])
# print(df_re[0].tolist()[:3])
# for sample in samples3:
# 	if sample[:-len('_bismark')] in df_re[0].tolist():
# 		sample_new = df_re.loc[df_re[0]==sample[:-len('_bismark')], 1].values[0]
# 		src = os.path.join(path3, sample)
# 		dst = os.path.join(path3, sample_new)
# 		# print(dst)
# 		os.rename(src, dst)

# for sample, sample_new in zip(samples1, samples1_new):
# 	if sample_new in df_m1.Sample.tolist():
# 		src = os.path.join(path1, sample)
# 		dst = os.path.join(path1, sample_new)
# 		os.rename(src, dst)

# metadata
# f_meta = './data/metadata/metadata_MB_v1_updated.tsv'
# df_meta = pd.read_table(f_meta, index_col='Sample')
# print(df_meta.head())


# allc 

# check if folders are empty
# allc_path = './data/allc/MB_v1'

# allcs = os.listdir(allc_path)
# i = 0
# for allc in allcs:
# 	if os.listdir(os.path.join(allc_path, allc)):
# 		i+=1
# print(i)
# print(len(allcs))


# check if allc matches with metadata

# rename files
path1 = './data/allc/MB_v1'
path2 = './data/allc/MB_EA'
path3 = './data/allc/MB_EB'

samples1 = os.listdir(path1)
samples2 = os.listdir(path2)
samples3 = os.listdir(path3)

print(len(samples1))
print(len(samples2))
print(len(samples3))



# match metadata
df_m1 = pd.read_table('./data/metadata/metadata_MB_v1_updated.tsv')
df_m2 = pd.read_table('./data/metadata/metadata_MB_EA_updated.tsv')
df_m3 = pd.read_table('./data/metadata/metadata_MB_EB_updated.tsv')

print(df_m1.shape)
print(df_m2.shape)
print(df_m3.shape)


print(len(set(samples1) & set(df_m1.Sample.tolist())))
print(len(set(samples2) & set(df_m2.Sample.tolist())))
print(len(set(samples3) & set(df_m3.Sample.tolist())))



# rename allc files in renamed allc folders 
path1 = './data/allc/MB_v1'
path2 = './data/allc/MB_EA'
path3 = './data/allc/MB_EB'

# set allc file names
def rename_allc(path1):
	samples1 = os.listdir(path1)
	for sample in samples1:
		pth = os.path.join(path1, sample)
		files = os.listdir(pth)
		gz_files = [file for file in files if file.endswith('.gz')]
		tbi_files = [file for file in files if file.endswith('.tbi')]
		for gz_file in gz_files:
			chrom=gz_file.split('_')[-1][:-len('.tsv.gz')]
			src = os.path.join(pth, gz_file)
			gz_file_new = 'allc_'+sample+'_'+chrom+'.tsv.gz'
			dst = os.path.join(pth, gz_file_new) 
			# print(src)
			# print(dst)
			os.rename(src, dst)

		for tbi_file in tbi_files:
			chrom=tbi_file.split('_')[-1][:-len('.tsv.gz.tbi')]
			src = os.path.join(pth, tbi_file)
			tbi_file_new = 'allc_'+sample+'_'+chrom+'.tsv.gz.tbi'
			dst = os.path.join(pth, tbi_file_new) 
			# print(src)
			# print(dst)
			os.rename(src, dst)
	return 0


# rename_allc(path1)
# check if path2 and path3 have the same naming convention with path1 
rename_allc(path2)
rename_allc(path3)

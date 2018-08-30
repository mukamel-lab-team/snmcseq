#!/usr/bin/env python3
import os
import glob
import pandas as pd
import numpy as np


def split_by_stride(np_array, stride):
	"""
	args: a np_array

	return: a list of np_arrays
	"""
	stride = int(stride)
	array_list = []
	start = 0
	end = start + stride 

	while end < len(np_array):
		array_list.append(np_array[start:end])			 
		start += stride
		end += stride

	if start < len(np_array):
		residue = np_array[start:len(np_array)]
		array_list.append(residue)

	return array_list


def get_sample_from_dir(local_dir):
	"""
	remove possible trailing "_bismark" from directory name
	"""
	dir_basename = os.path.basename(local_dir) 
	if dir_basename.endswith('_bismark'):
		sample_name = dir_basename[:-len('_bismark')]		
	else:
		sample_name = dir_basename

	return sample_name	

def gen_tsne_input(local_dirs, 
		output_filename='./tsne/tsne_input.tsv', 
		pooling_num=None):
	"""
	Generate the input file required by the downstream tSNE analysis
	from methylation profiles of samples   

	args: a list of local_dir paths containing methylation profiles 
		(each local_dir stores the methylome of a sample cell)

	output: a tsv matrix containing _mc and _c columns for each sample
		(CH methylation)
	"""

	df_combined = pd.DataFrame()
	# loop over different samples (local_dirs)
	tc = []
	for i, local_dir in enumerate(sorted(local_dirs)):
		sample_name = get_sample_from_dir(local_dir)	
		print("Processing " + sample_name + "...")
		# loop over different bin files within a sample
		sample_mch = np.array([], dtype=int) 
		sample_ch = np.array([], dtype=int)
		sample_mcg = np.array([], dtype=int) 
		sample_cg = np.array([], dtype=int)
		# i==0 dealing with chr and bin
		if i == 0:
			chrom = np.array([], dtype=np.unicode)
			bin_pos = np.array([], dtype=int)
		else:
			chrom_test = np.array([], dtype=np.unicode)
			bin_pos_test = np.array([], dtype=int)

		# loop over different chromomsomes
		# chromosome order is very important!!!!
		for fname in sorted(glob.glob(os.path.join(local_dir, '*_10000_*.tsv'))):
			# 0/1/2/3/4/5/: chr/bin/mcg/cg/mch/ch	
			# process 1 chromosome
			df = pd.read_table(fname, header=None)

			# i==0 dealing with chr and bin
			if i == 0:
				chr_chrom = df[0].values
				chr_bin_pos = df[1].values
				if pooling_num:
					chr_chrom = [array[0] 
							for array in 
							split_by_stride(chr_chrom, pooling_num)]	
					chr_bin_pos = [array[0] 
							for array in 
							split_by_stride(chr_bin_pos, pooling_num)]	

				chrom = np.concatenate((chrom, chr_chrom))
				bin_pos = np.concatenate((bin_pos, chr_bin_pos))

			else:
				chr_chrom = df[0].values
				chr_bin_pos = df[1].values
				if pooling_num:
					chr_chrom = [array[0] 
							for array in 
							split_by_stride(chr_chrom, pooling_num)]	
					chr_bin_pos = [array[0] 
							for array in 
							split_by_stride(chr_bin_pos, pooling_num)]	

				chrom_test = np.concatenate((chrom_test, chr_chrom))
				bin_pos_test = np.concatenate((bin_pos_test, chr_bin_pos))

			# dealing with mCH and CH 
			sample_chr_mch = df[4].values 
			sample_chr_ch = df[5].values 
			sample_chr_mcg = df[2].values 
			sample_chr_cg = df[3].values 
			# sum over neighboring bins at a certain "pooling rate"
			if pooling_num:
				sample_chr_mch = [np.sum(array) 
							for array in 
							split_by_stride(sample_chr_mch, pooling_num)]	
				sample_chr_ch = [np.sum(array) 
							for array in 
							split_by_stride(sample_chr_ch, pooling_num)]	
				sample_chr_mcg = [np.sum(array) 
							for array in 
							split_by_stride(sample_chr_mcg, pooling_num)]	
				sample_chr_cg = [np.sum(array) 
							for array in 
							split_by_stride(sample_chr_cg, pooling_num)]	

			sample_mch = np.concatenate((sample_mch, sample_chr_mch)) 
			sample_ch = np.concatenate((sample_ch, sample_chr_ch))	

			sample_mcg = np.concatenate((sample_mcg, sample_chr_mcg)) 
			sample_cg = np.concatenate((sample_cg, sample_chr_cg))	

		if i != 0:
			# check the order of chromosomes and pos 
			assert chrom.tolist() == chrom_test.tolist()
			assert bin_pos.tolist() == bin_pos_test.tolist()
		# assert len(chrom) == len(bin_pos) == len(sample_ch) == len(sample_mch)
		if not (len(chrom) == len(bin_pos) 
			== len(sample_ch) == len(sample_mch)
			== len(sample_cg) == len(sample_mcg)):
			print("Corrupted sample: %s, removed from dataset" % sample_name)
			continue	
		# ignore empty samples
		if len(sample_mch) == 0: 
			print("Empty sample: %s" % sample_name)
			continue

		# output
		if i == 0:
			col_names = ['chr', 'bin']
			df_combined['chr'] = chrom
			df_combined['bin'] = bin_pos 

			tc = sample_ch + sample_cg
		else:
			tc += (sample_ch+sample_cg)

		# col_names.append(sample_name+'_mch')
		# col_names.append(sample_name+'_ch')
		# col_names.append(sample_name+'_mcg')
		# col_names.append(sample_name+'_cg')
		# df_combined[sample_name+'_mch'] = sample_mch	
		# df_combined[sample_name+'_ch'] = sample_ch 
		# df_combined[sample_name+'_mcg'] = sample_mcg	
		# df_combined[sample_name+'_cg'] = sample_cg 

		# col_names.append(sample_name+'_tc')
		# df_combined[sample_name+'_tc'] = sample_ch + sample_cg	


	mean_tc = tc/len(local_dirs)
	col_names.append('mean_tc')
	df_combined['mean_tc'] = mean_tc 

	# output format:
	# columns: chr, bin, sample1_mch, sample1_ch, sample1_mcg, sample1_cg ...
	df_combined.reindex_axis(col_names, axis=1)
	df_combined.to_csv(output_filename, sep='\t', index=False)
	print(df_combined.head())
	print("All done!")
	return 0


if __name__ == '__main__':

	DIR = '/cndd/Public_Datasets/single_cell_methylome/binc/human_combined'
	local_dirs = glob.glob(os.path.join(DIR, '*'))
	# output_filename = '/cndd/Public_Datasets/single_cell_methylome/binc/human_combined_100kb_bins.tsv' 
	output_filename = '/cndd/fangming/side_project/snmcseq_dev/binc/human_combined_100kb_bins_mean_tc.tsv' 
	# DIR = '/cndd/Public_Datasets/single_cell_methylome/binc/human'
	# local_dirs = glob.glob(os.path.join(DIR, 'Pool_*'))
	# output_filename = './tsne/human_combined_100000.tsv' 


	gen_tsne_input(local_dirs, output_filename=output_filename, pooling_num=10)


# df_combined = pd.DataFrame()

# for fname in fnames:
# 	sample_name = os.path.basename(fname).split('.')[0]
# 	df = pd.read_table(fname, header=None)
# 	df = df[[4,5]].head(100)

# 	df_combined[sample_name+'_mc'] = df[4].values	
# 	df_combined[sample_name+'_c'] = df[5].values 

# # output format:
# # columns: sample1_mc, sample1_c, sample2_mc, sample2_c, ...
# df_combined.to_csv('test_tsne.tsv', sep='\t', index=False)

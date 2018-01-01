#!/usr/bin/env python3

import argparse 
import pandas as pd

def create_parser():
	"""
	"""
	parser = argparse.ArgumentParser(description="Outputs the TSNE coordinates for a set of input cells. Input is the methylated and total base calls for each cell across bins/genes in tsv format. "+
	                                 "DETAILS: Load the data, filter out any samples not in the QC list, impute mCH at bins/gene where coverage is poor, normalize the data, "+
	                                 "run PCA (50 components) and TSNE (2 components) and save the output to file.")

	parser.add_argument("-i", "--input", help="input file name. Each column is a cluster label",
						required=True)
	parser.add_argument("-m", "--mdata", help="path to metdata file",
						required=True)
	parser.add_argument("-o", "--output", help="output file name", 
						required=True)
	return parser		


def get_clusters(input_fname):
	"""
	get clusters from the output of backspin.py 

	args: input file name format: each line is a cluster

	return: a dict
	""" 

	with open(input_fname, 'r') as input_file:
		line = input_file.readline()	
		label = 1
		cluster_dict = dict()
		while line:
			samples = line.strip('\n').split(',')
			cluster_dict['cluster_'+str(label)] = [sample.strip('_mcc') for sample in samples]

			label += 1
			line = input_file.readline()

	# print(len(cluster_dict))
	return cluster_dict


def dict_to_dataframe(cluster_dict, sample_list):
	"""
	transform a cluster dictionary to a pandas dataframe
	"""


	# iterate over all clusters
	dict_list = []
	for cluster_label, samples in cluster_dict.items():
		# iterate over all samples
		for sample in samples:
			sample_dict = {'cluster_ID': cluster_label, 'sample': sample}
			dict_list.append(sample_dict)

	df = pd.DataFrame(dict_list)

	# for sample in sample_list:
	# 	if sample not in df['sample'].tolist():
	# 		dict_list.append({'cluster_ID': 'NA', 'sample': sample})
	# df = pd.DataFrame(dict_list)	

	df = df[['sample', 'cluster_ID']]
	df.sort_values('sample', inplace=True)
	df = df.reset_index(drop=True)
	assert df['sample'].tolist() == sample_list
	return df


def get_sample_list(metadata_fname):
	"""
	"""
	df_meta = pd.read_table(metadata_fname, header=0)
	return sorted(df_meta['Sample'].tolist())


if __name__ == '__main__':
	parser = create_parser()
	args = parser.parse_args()
	cluster_dict = get_clusters(args.input)
	sample_list = get_sample_list(args.mdata)
	cluster_df = dict_to_dataframe(cluster_dict, sample_list)
	cluster_df.to_csv(args.output, sep='\t', header=True, index=False)

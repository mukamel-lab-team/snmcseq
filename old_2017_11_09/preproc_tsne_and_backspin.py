#!/usr/bin/env python3

import numpy as np
import pandas as pd
from collections import OrderedDict
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import ipdb
import random
import math
import numpy as np
import sys, getopt
import argparse
import warnings

# import mypy

def get_sample_names(df):
	"""
	return: a list of sample names
	"""
	col_names = df.columns.tolist()

	samples = []
	for col_name in col_names:
		if col_name.endswith('_c'):
			col_name = col_name[:-len('_c')] 
		elif col_name.endswith('_mc'):
			col_name = col_name[:-len('_mc')]
		else:
			raise 	
		samples.append(col_name)

	sample_names = list(set(samples))
	return sample_names


warnings.filterwarnings("ignore")

###################################
# Parse command line arguments
###################################
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                 description="Outputs the TSNE coordinates for a set of input cells. Input is the methylated and total base calls for each cell across bins/genes in tsv format. "+
                                 "DETAILS: Load the data, filter out any samples not in the QC list, impute mCH at bins/gene where coverage is poor, normalize the data, "+
                                 "run PCA (50 components) and TSNE (2 components) and save the output to file.")

parser.add_argument("-i", "--input", help="input file name. Must have columns named as sample1_mc, sample1_c, sample2_mc, etc.", 
                    required=True)
# added 10/31/2017 ---
parser.add_argument("-ex", "--ex_cols", help="number of the first few columns to remove", 
                    type=int, required=True)
parser.add_argument("-cx", "--context", help="CG or CH", 
                    required=True)
# end of the adding --- 

parser.add_argument("-o", "--output", help="output file name", required=True)
# parser.add_argument("-s", "--species", help="mouse or human", default="mouse")
parser.add_argument("-n", "--normalize", help="normalize the data before running PCA and TSNE", action='store_true')
# parser.add_argument("-p", "--perplexity", type=int, help="TSNE perplexity", default=25)
# parser.add_argument("-d", "--seed", type=int, help="TSNE seed", default=1)
parser.add_argument("-b", "--base_call_cutoff", type=int, help="minimum base calls for a bin to not be imputed.", default=100)
parser.add_argument("-m", "--mdata", help="path to metdata file")
args = parser.parse_args()

# species = args.species
normalize = args.normalize
# perplexity = args.perplexity
infile = args.input
outfile = args.output
base_call_cutoff = args.base_call_cutoff
# seed = args.seed
mdata = args.mdata
ex_cols = args.ex_cols
context = args.context


##################################
# Load data
##################################
print("Loading data.")



# Load input and metadata
df_gene_level_mCH = pd.read_csv(infile, sep="\t")
metadata = pd.read_csv(mdata, sep="\t")

# df = df_gene_level_mCH
df = df_gene_level_mCH.iloc[:, ex_cols:] ##### remove chrom and bin position

# samples = df.samples.tolist()
samples = get_sample_names(df)

print("Computing m%s levels." % context)

# Keep only bins that have sufficient coverage in at least 99.5% of all cells
df = df.loc[(df.filter(regex='_c$') > base_call_cutoff).sum(axis=1) >= .995*len(samples)]
print("Matrix size after pruning... "+ str(df.shape))
df_mc = df[[x+'_mc' for x in samples]]
df_c = df[[x+'_c' for x in samples]]
df_c[df_c < base_call_cutoff] = np.nan
df_c.columns = [x+'_mcc' for x in samples]
df_mc.columns = [x+'_mcc' for x in samples]
df = df_mc/df_c

print('Imputing data.')
# Impute missing values
df = df.loc[df.count(axis=1) > 0]
df.reset_index(inplace=True, drop=True)
means = df.mean(axis=1)
fill_value = pd.DataFrame({col: means for col in df.columns})
df.fillna(fill_value, inplace=True)

if normalize:
    print('Normalizing.')
   	# check if samples in df match samples in metadata
    # fangming 09/09/2017	 
    samples_filtered = []
    # extract a sample only if it's in the metadata
    for sample in samples:
        if sample in metadata['Sample'].values:
            samples_filtered.append(sample) 
    # truncate dfs fangming 09/09/2017
    df_columns_filtered = [item+'_mcc' for item in samples_filtered]
    # print(df.shape)
    df = df[df_columns_filtered]
    # print(df.shape)

    for i,row in metadata.iterrows():
        samp = row['Sample']
        if samp+'_mcc' in df.columns:  # fangming edit 09/04/2017
            if context == 'CH':
                df[samp+'_mcc'] = (df[samp+'_mcc'] / (row['mCH/CH']+.01))
            elif context == 'CG':
                df[samp+'_mcc'] = (df[samp+'_mcc'] / (row['mCG/CG']+.01))
            else:
                raise ValueError('Wrong context: %s' % context)

df.to_csv(outfile, sep='\t', header=True, index=False)

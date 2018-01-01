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
import time

ti = time.time()
warnings.filterwarnings("ignore")

###################################
# Parse command line arguments
###################################
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                 description="Outputs the TSNE coordinates for a set of input cells. Input is the methylated and total base calls for each cell across bins/genes in tsv format. "+
                                 "DETAILS: Load the data, filter out any samples not in the QC list, impute mCH at bins/gene where coverage is poor, normalize the data, "+
                                 "run PCA (50 components) and TSNE (2 components) and save the output to file.")

parser.add_argument("-i", "--input", help="input file name. Columns should be named as sample1_mcc, sample2_mcc, etc.", 
                    required=True)
parser.add_argument("-o", "--output", help="output file name", required=True)
parser.add_argument("-p", "--perplexity", type=int, help="TSNE perplexity", default=25)
parser.add_argument("-d", "--seed", type=int, help="TSNE seed", default=1)
args = parser.parse_args()

perplexity = args.perplexity
infile = args.input
outfile = args.output
seed = args.seed

#######################
# LOADING DATA (preprocessed)
#######################
print('Loading data...')
df = pd.read_table(infile) 
df = df.filter(regex='_mcc$')

print('Data matrix shape:')
print(df.shape)

#######################
# RUNNING PCA
#######################
print("Running PCA.")
pca = PCA(n_components=50)
sklearn_transf = pca.fit_transform(df.T)
# print(pca.explained_variance_ratio_)
sklearn_transf_PCA = sklearn_transf



#######################
# RUNNING TSNE
#######################
print("Running TSNE.")
num_components = 2
tsne = TSNE(n_components=num_components, init='pca', random_state=seed, perplexity=perplexity, verbose=3)
sklearn_transf = tsne.fit_transform(sklearn_transf_PCA)

print("Saving output to file.")
df_tsne = pd.DataFrame(sklearn_transf, columns=['tsne_x','tsne_y'])
df_tsne['sample'] = [sample[:-len('_mcc')] for sample in df.columns.tolist()]
df_tsne.to_csv(outfile, sep="\t", na_rep='NA', header=True, index=False)


tf = time.time()
print("Running time: %.2f seconds." % (tf - ti))
# plt.plot(sklearn_transf[:,0], sklearn_transf[:,1], '.r')


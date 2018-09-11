#!/usr/bin/env python3

import matplotlib
matplotlib.use('Agg')

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import itertools
import time
import os
import logging
import matplotlib.gridspec as gridspec

from sklearn.decomposition import PCA
from sklearn.manifold import TSNE 

from __init__ import *
from snmcseq_utils import create_logger

BIN_COV_THRESHOLD = 0.8 # a bin will be included if it has coverage in over #(frac) of cells.

def default_paras(context):
	"""
	"""
	if context == 'CH':
		base_call_cutoff = 250
		fraction_included = 0.8
		high_mean = 4000
	elif context == 'CA':
		base_call_cutoff = 100
		fraction_included = 0.8
		high_mean = 1400
	elif context == 'CG':
		base_call_cutoff = 10
		fraction_included = 0.8
		high_mean = 300
	else:	
		raise ValueError('Incorrect context, choose from {}'.format(CONTEXTS))

	return base_call_cutoff, fraction_included, high_mean

def normalize_by_global(df_mcc, df_meta, context):
	"""
	normalize a mcc dataframe by global mcc level
	"""
	logging.info('Normalize by global m%s ...' % context)
	df_nmcc = pd.DataFrame()

	for idx, row in df_meta.iterrows():
		try:
			samp = idx 
			if context == 'CH':
				df_nmcc[samp] = (df_mcc[samp] / (row['mCH/CH']+.01))
			elif context == 'CG':
				df_nmcc[samp] = (df_mcc[samp] / (row['mCG/CG']+.01))
			elif context == 'CA':
				df_nmcc[samp] = (df_mcc[samp] / (row['mCA/CA']+.01))
			else:
				raise ValueError('Incorrect context, choose from {}'.format(CONTEXTS))
		except:
			logging.warning("Warning: {} in mapping summary not found in df_mcc".format(samp))

	return df_nmcc


def preproc_summary(num_kept, num_orig, frac_kept, lower_limit, num_high_removed, output_summary=None):
	# summary
	summary = ( 
	'''***** Preprocessing summary *****
	Number of bins kept: {}/{} ({})
	Lower limit of coverage: {}
	Number of removed high coverage bins: {}
	***** End of preprocessing summary *****
	''').format(
	    num_kept, num_orig, frac_kept, 
	    lower_limit, num_high_removed)

	logging.info(summary)
	if output_summary:
		with open(output_summary, 'w') as file:
		    file.write(summary)
		logging.info("Saved QC summary to {}".format(output_summary))
	return summary

def preproc_plot(ens, df_c, condition, fcells, fcells_after, fbins, fbins_after, context, output_plot=None):
	"""
	"""
	# plot
	if context in ['CH', 'CA']:
		zoom1 = 5000
		zoom2 = 1500
	elif context == 'CG':
		zoom1 = 500
		zoom2 = 150
	else:
		raise ValueError('Incorrect context, choose from {}'.format(CONTEXTS))

	# evaluate the filtering process (what're kept)
	baseline = df_c.mean(axis=1)
	exp = condition.values*baseline.values

	fig = plt.figure(figsize=(12, 15))
	gs = gridspec.GridSpec(5, 2)
	ax = fig.add_subplot(gs[0, :])
	ax.set_title('binc_m{}_QC_{}'.format(context, ens))
	ax.plot(baseline.values, linewidth=0.5, label='All bins')
	ax.plot(exp, 'o', markersize=1, label='Included bins')
	ax.legend()
	ax.set_ylabel('Total cytosine coverage')

	ax = fig.add_subplot(gs[1, :])
	ax.plot(baseline.values, linewidth=0.5, label='All bins')
	ax.plot(exp, 'o', markersize=1, label='Included bins')
	ax.legend()
	ax.set_ylim([0, zoom1])
	ax.set_xlabel('100kb bins across mm10 genome')
	ax.set_ylabel('Total cytosine coverage')

	ax = fig.add_subplot(gs[2, :])
	ax.plot(baseline.values, linewidth=0.5, label='All bins')
	ax.plot(exp, 'o', markersize=1, label='Included bins')
	ax.legend()
	ax.set_ylim([0, zoom2])
	ax.set_xlabel('100kb bins across mm10 genome')
	ax.set_ylabel('Total cytosine coverage')

	ax = fig.add_subplot(gs[3:, 0])
	ax.plot(fbins.sort_values().values, label='Before')
	ax.plot(fbins_after.sort_values().values, label='After')
	ax.legend()
	ax.set_title('Bins QC')
	ax.set_xlabel('bins')
	ax.set_ylabel('Fraction of cells passing cutoff coverage for that bin')

	ax = fig.add_subplot(gs[3:, 1])
	ax.plot(fcells.sort_values().values, label='Before')
	ax.plot(fcells_after.sort_values().values, label='After')
	ax.legend()
	ax.set_title('Cells QC')
	ax.set_xlabel('cells')
	ax.set_ylabel('Fraction of bins passing cutoff coverage for that cell')

	fig.tight_layout()
	if output_plot:
		fig.savefig(output_plot)
		logging.info("Saved QC plot to {}".format(output_plot))
		plt.close('all')
	else:
		plt.show()

	return

def calculate_nmcc(df_mc, df_c_nan, df_meta, condition, context, output_nmcc=None):

	df_mcc = df_mc[condition]/df_c_nan[condition]
	# normalization (normalized by global mcc)
	logging.info('Normalization...')
	df_mcc = normalize_by_global(df_mcc, df_meta, context)

	# imputation (missing value -> mean value of all cells)
	logging.info('Imputing data...')
	means = df_mcc.mean(axis=1)
	fill_value = pd.DataFrame({col: means for col in df_mcc.columns})
	df_mcc.fillna(fill_value, inplace=True)

	# add "_mcc" suffix
	df_mcc.columns = df_mcc.columns.values + '_mcc'
	logging.info('Done preprocess bins...\nOutput shape: {}'.format(df_mcc.shape))

	# save to output
	if output_nmcc:
		df_mcc.to_csv(output_nmcc, sep='\t', na_rep='NA', header=True, index=True)
		logging.info('Saved normalized mcc file to: {}'.format(output_nmcc))

	return df_mcc


def preproc_cellcovr(fcells, fcells_after, output_cellcovr=None):
	# files (cell conf level)
	cells_covr = pd.merge(fcells.to_frame(), fcells_after.to_frame(), left_index=True, right_index=True)
	cells_covr.columns=['covr_before', 'covr_after']
	cells_covr = cells_covr.rename_axis('sample')
	if output_cellcovr:
		cells_covr.to_csv(output_cellcovr, sep='\t', na_rep='NA', header=True, index=True)
		logging.info("Saved cells coverage info to {}".format(output_cellcovr))
	return cells_covr


def preproc_bins(ens, df_meta=None, df=None, bin_size=BIN_SIZE_FEATURE, context='CH', base_call_cutoff=None, fraction_included=None, high_mean=None, to_file=False):
	"""
	df is df_mc_c
	df and df_meta could also be got from pipeline directly

	recommended parameters:
	mCH: 250, 0.8, 4000
	mCA: 100, 0.8, 1400 (1600)	
	mCG: 10, 0.8, 300 (140)
	"""
	# begin preproc
	logging.info("Begin preprocessing binc data...{}_{}_{}".format(ens, bin_size, context))

	# set up paths
	ens_path = os.path.join(PATH_ENSEMBLES, ens)
	output_summary = os.path.join(ens_path, 'binc_m{}_{}_QC_summary_{}.tsv'.format(context, bin_size, ens))
	output_plot = os.path.join(ens_path, 'plots/binc_m{}_{}_QC_plot_{}.pdf'.format(context, bin_size, ens)) 
	output_cellcovr = os.path.join(ens_path, 'binc_m{}_{}_QC_cells_covr_{}.tsv'.format(context, bin_size, ens))

	if to_file:
		output_nmcc = os.path.join(ens_path, 'binc/binc_m{}_{}_nmcc_{}.tsv'.format(context, bin_size, ens))
	else:
		output_nmcc = None 


	if not os.path.exists(os.path.join(ens_path, 'plots')):
	    os.makedirs(os.path.join(ens_path, 'plots'))

	# data loading 
	if not isinstance(df_meta, pd.DataFrame):
		meta_file = os.path.join(ens_path, 'mapping_summary_{}.tsv'.format(ens))
		# metadata
		df_meta = pd.read_table(meta_file, index_col='Sample')
	if not isinstance(df, pd.DataFrame):
		binc_file = os.path.join(ens_path, 'binc/binc_m{}_{}_{}.tsv.bgz'.format(context, bin_size, ens))
		df = pd.read_table(binc_file, index_col=['chr', 'bin'], 
			compression='gzip', dtype={'chr': object})

	# set up parameters
	if not (base_call_cutoff and fraction_included and high_mean):
		base_call_cutoff, fraction_included, high_mean = default_paras(context)


	# get df_mc and df_c   
	df_mc = df.filter(regex='_mc$')
	df_c = df.filter(regex='_c$')
	df_c.columns = [col[:-len('_c')] for col in df_c.columns] 
	df_mc.columns = [col[:-len('_mc')] for col in df_mc.columns] 


	# low coverage condition
	df_c_nan = df_c.copy()
	df_c_nan[df_c < base_call_cutoff] = np.nan
	nonnulls = 1 - df_c_nan.isnull()
	fbins = nonnulls.sum(axis=1)/nonnulls.shape[1]
	fcells = nonnulls.sum(axis=0)/nonnulls.shape[0]
	# rank bins by fbins and get top fraction for each chromosome 
	dfs = []
	for chrom, sr_sub in fbins.groupby('chr'):
	    num_bins = sr_sub.shape[0]
	    candicates = sr_sub.nlargest(int(fraction_included*num_bins))
	    # Fangming 09/11/2018 to remove nan values
	    candidates = candidates[candidates>BIN_COV_THRESHOLD]
	    dfs.append(candidates)
	dfs = pd.concat(dfs)    
	condition_low = fbins.index.isin(dfs.index) # just for the purpose of re-ordering index 

	# high coverage condition
	condition_high = df_c.mean(axis=1) < high_mean 
	# combined condition
	condition = condition_high & condition_low               

	# status after preprocessing
	lower_limit = fbins[condition_low].nsmallest(10).mean()
	num_high_removed = np.sum(~condition_high)
	num_kept = fbins[condition].shape[0]
	num_orig = fbins.shape[0]
	frac_kept = num_kept/num_orig
	fbins_after = nonnulls[condition].sum(axis=1)/nonnulls[condition].shape[1]
	fcells_after = nonnulls[condition].sum(axis=0)/nonnulls[condition].shape[0]

	# outputs
	# summary
	summary = preproc_summary(num_kept, num_orig, frac_kept, lower_limit, num_high_removed, output_summary=output_summary)
	# cellcovr
	cellcovr = preproc_cellcovr(fcells, fcells_after, output_cellcovr=output_cellcovr)	
	# plot
	preproc_plot(ens, df_c, condition, fcells, fcells_after, fbins, fbins_after, context, output_plot=output_plot)
	# nmcc
	df_mcc = calculate_nmcc(df_mc, df_c_nan, df_meta, condition, context, output_nmcc=output_nmcc)

	return df_mcc


def preproc_bins_combine_contexts(ens, dfs=None, bin_size=BIN_SIZE_FEATURE, contexts=['CH', 'CG'], to_file=False):
	"""Concatenate and normalize (by standard deviation) different marks of bins
	given df_nmcc files, concatenate and normalize them

	"""
	assert len(contexts) == 2

	# begin preproc
	logging.info("Begin preprocessing binc data...{}_{}_{}".format(ens, bin_size, contexts))

	# set up paths
	ens_path = os.path.join(PATH_ENSEMBLES, ens)
	if to_file:
		output_nmcc = os.path.join(ens_path, 'binc/binc_m{}m{}_{}_nmcc_{}.tsv'.format(contexts[0], contexts[1], bin_size, ens))
	else:
		output_nmcc = None 

	# data loading 
	if (not dfs) or (len(dfs) != 2):
		dfs = []
		for context in contexts:
			binc_file = os.path.join(ens_path, 'binc/binc_m{}_{}_nmcc_{}.tsv'.format(context, bin_size, ens))
			df = pd.read_table(binc_file, index_col=['chr', 'bin'], 
				# compression='gzip', 
				dtype={'chr': object})
			dfs.append(df)

	# get an overall standard deviation
	stds = [df.stack().std() for df in dfs]
	logging.info("Standard deviation before normalization: {}-{}".format(contexts, stds))

	# normalize by over all standard deviation
	dfs_n = [df/std for df, std in zip(dfs, stds)]
	stds_n = [df_n.stack().std() for df_n in dfs_n]
	logging.info("Standard deviation after normalization: {}-{}".format(contexts, stds_n))

	# concatenate
	df_cat = pd.concat(dfs_n)
	logging.info("Input shapes: {}_{}\nOutput shape: {}".format(dfs[0].shape, dfs[1].shape, df_cat.shape))

	# output
	if to_file:
		df_cat.to_csv(output_nmcc, sep='\t', na_rep='NA', header=True, index=True)
		logging.info("Saved combined feature matrix to {}".format(output_nmcc))

	return df_cat


if __name__ == '__main__':
	log = create_logger()


	# configs 

	# ens = 'Ens3'
	# context = 'CH'
	ens = 'Ens33'
	for context in CONTEXTS:
		preproc_bins(ens, bin_size=BIN_SIZE_FEATURE, context=context, to_file=True)

	ens = 'Ens33'
	for contexts in COMBINED_CONTEXTS_LIST:
		preproc_bins_combine_contexts(ens, dfs=None, bin_size=BIN_SIZE_FEATURE, contexts=contexts, to_file=True)

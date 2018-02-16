#!/usr/bin/env python3
# coding: utf-8

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
import os
import glob
import logging
import time

from scipy.stats import spearmanr
from scipy.stats import zscore 

from __init__ import *
from snmcseq_utils import create_logger 
from snmcseq_utils import plot_tsne_labels
from snmcseq_utils import create_logger 

OUTPUT_DIR = '/cndd/fangming/snmcseq_dev/data/cluster/test_clusters/MB_v1'

# reorder columns
def reorder_cols(cols):
    cols_ordered = sorted([int(col.split('_')[1]) for col in cols])
    cols_ordered = ['cluster_'+str(col)+'_mcc' for col in cols_ordered]
    return cols_ordered

def auto_annotate_worker(cluster_f, input_f, 
    output_f, output_heatmap, 
    reject_threshold = 0.9,
    input_fgt=os.path.join(PATH_REFERENCES, 'Mouse_published/binc_mCH_100000_clusterwise_mcc_mouse_published.tsv'), 
    input_gt=os.path.join(PATH_REFERENCES, 'Mouse_published/summary_mouse_published.tsv')):
    """
    """
    ti = time.time()
    logging.info('Auto annotating {}'.format(cluster_f))

    # load data
    df_input = pd.read_table(input_f, 
        index_col=['chr', 'bin'], dtype={'chr': object}, compression='gzip')
    df_cluster = pd.read_table(cluster_f, index_col='sample')


    # cluster mc_c
    df_c = df_input.filter(regex='_c$')
    df_mc = df_input.filter(regex='_mc$')

    df_mc_c = pd.DataFrame() 
    for label, df_sub in df_cluster.groupby('cluster_ID'):
        samples = df_sub.index.values
        df_mc_c[label+'_mc'] = df_mc[samples+'_mc'].sum(axis=1)
        df_mc_c[label+'_c'] = df_c[samples+'_c'].sum(axis=1)

    # cluster mcc
    cov = (df_mc_c.filter(regex='_c$').apply(np.log10)>3)*1
    gene_cov = cov.apply(np.all, axis=1)
    df_filtered = df_mc_c.loc[gene_cov, :]

    df_mc = df_filtered.filter(regex='_mc$')
    df_c = df_filtered.filter(regex='_c$')
    df_mc.columns = [col[:-len('_mc')] for col in df_mc.columns]
    df_c.columns = [col[:-len('_c')] for col in df_c.columns]
    df_mcc = df_mc/df_c
    df_mcc.columns = [col+'_mcc' for col in df_mcc.columns]
    # logging.info(df_mcc.shape)

    # compare with ground truth
    df = pd.read_table(input_fgt, 
        index_col=['chr', 'bin'], dtype={'chr': object})
    df2 = df_mcc 
    df_gt = pd.read_table(input_gt, index_col='Sample')['Neuron type'].to_frame()

    # combine gt clusters and clusters
    df_cmb = pd.merge(df, df2, left_index=True, right_index=True)
    # logging.info(df_cmb.shape)

    # reorder columns
    celltype_cols = df.filter(regex='_mcc$').columns.values 
    cluster_cols = reorder_cols(df2.filter(regex='_mcc$').columns.values) 
    celltypes = [col[:-len('_mcc')] for col in celltype_cols]
    clusters = [col[:-len('_mcc')] for col in cluster_cols]

    # do spearman correlation cell_type * clusters
    rho, pval = spearmanr(df_cmb[celltype_cols], df_cmb[cluster_cols])
    corr = rho[:len(celltype_cols), len(celltype_cols):]
    df_corr = pd.DataFrame(corr, columns=clusters, index=celltypes)
    df_corr_zscore = df_corr.apply(zscore, axis=0)

    df_annot = pd.DataFrame(columns=df_corr.columns, index=df_corr.index)
    celltype_res = []
    for col, idx in df_corr.idxmax().iteritems():
        if df_corr.loc[idx, col] > reject_threshold:
            # print(col + ' --> ' + idx)
            celltype_res.append({'cluster_ID': col,
                                 'cluster_annotation': idx})
            df_annot.loc[idx, col] = '*'
        else:
            # print(col + ' cell type not matched!')
            celltype_res.append({'cluster_ID': col,
                                 'cluster_annotation': np.nan})
            
    celltype_res = pd.DataFrame(celltype_res)
    celltype_res.to_csv(output_f, sep='\t', header=True, index=False, na_rep='NA')
    logging.info('Saved autotation file to {}'.format(output_f))

    df_annot = df_annot.fillna('')
       
    # plot 
    fig, axs = plt.subplots(2, 1, figsize=(8, 10))
    ax = axs[0]
    sns.heatmap(df_corr, 
            ax=ax,
            cmap='viridis',
            annot = df_annot, fmt='',
            # linewidth = 0.01,
            xticklabels = False,
            cbar_kws={'label':'Spearmanr correlation'})
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
    ax.set_title('Correlations of 100kb bin mCH \n between known cell types and $de$ $novo$ found clusters')

    ax = axs[1]
    sns.heatmap(df_corr_zscore, 
            ax=ax,
            cmap='viridis',
            annot = df_annot, fmt='',
            # linewidth = 0.01,
            xticklabels = False,
            cbar_kws={'label':'Zscore across rows'})
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
    # ax.set_title('')
    fig.tight_layout()
    fig.savefig(output_heatmap)
    logging.info('Saved heatmap to {}'.format(output_heatmap))

    tf = time.time()
    logging.info('Done annotating {}, total time spent: {} second'.format(cluster_f, tf-ti))

def run_autoannotate_CEMBA(ens):
    """
    """
    ens_path = os.path.join(PATH_ENSEMBLES, ens)

    input_f = os.path.join(ens_path, 'binc/binc_mCH_100000_{}.tsv.bgz'.format(ens)) 
    cluster_files = sorted(glob.glob(os.path.join(ens_path, 'cluster/cluster_*_{}.tsv'.format(ens)))) 
    for cluster_f in cluster_files:
        output_heatmap = os.path.join(ens_path, 'plots/autoannot_heatmap_{}.pdf'.format(os.path.basename(cluster_f)[:-len('.tsv')])) 
        output_f = os.path.join(ens_path, 'cluster/{}.annot'.format(os.path.basename(cluster_f))) 
        auto_annotate_worker(cluster_f, input_f, output_f, output_heatmap) 
        plt.close('all')
    return

if __name__ == '__main__':

    log = create_logger()
    enss = ['Ens1', 'Ens2', 'Ens3', 'Ens4']
    for ens in enss:
        autoannotate_CEMBA(ens)

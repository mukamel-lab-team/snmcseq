#!/usr/bin/env python3
# coding: utf-8

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import glob

from scipy.stats import spearmanr
from scipy.stats import zscore 

from snmcseq_utils import plot_tsne_labels
from snmcseq_utils import create_logger 

OUTPUT_DIR = '/cndd/fangming/snmcseq_dev/data/cluster/test_clusters/MB_v1'
plt.ioff()

# reorder columns
def reorder_cols(cols):
    cols_ordered = sorted([int(col.split('_')[1]) for col in cols])
    cols_ordered = ['cluster_'+str(col)+'_mcc' for col in cols_ordered]
    return cols_ordered

def auto_annotate(cluster_fname, input_f, OUTPUT_DIR):
    """
    """
    print('Processing {}'.format(cluster_fname))
    # load data
    cluster_f = os.path.join(OUTPUT_DIR, cluster_fname) 
    #
    output_f = cluster_f + '.annot' 
    output_tsne = os.path.join(OUTPUT_DIR, 'figures/{}.tsne.pdf'.format(cluster_fname)) 
    output_heatmap = os.path.join(OUTPUT_DIR, 'figures/{}.heatmap.pdf'.format(cluster_fname)) 

    df_input = pd.read_table(input_f, index_col='id')
    df_cluster = pd.read_table(cluster_f)


    # cluster mc_c
    df_c = df_input.filter(regex='_c$')
    df_mc = df_input.filter(regex='_mc$')

    df_mc_c = pd.DataFrame() 
    for label, df_sub in df_cluster.groupby('cluster_ID'):
        samples = df_sub['sample'].values
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
    print(df_mcc.shape)

    # load ground truth
    input_fgt = './data/cluster/cluster_MB_v1/genebody_mCH_human_combined_cluster_MB_v1_mcc.tsv'
    input_tsne = './data/tsne/tsne_perp30_binc_mCH_human_combined_100000_summary_nmcc_v3.tsv'
    input_meta = './data/metadata/metadata_human_combined_updated.tsv'
    input_gt = './data/cluster/cluster_MB_v1/cluster_MB_v1.tsv'

    df = pd.read_table(input_fgt, index_col='id')
    df2 = df_mcc 
    df_meta = pd.read_table(input_meta, index_col='Sample')
    df_tsne = pd.read_table(input_tsne)
    df_gt = pd.read_table(input_gt, index_col='sample')['cell_type'].to_frame()

    # combine gt clusters and clusters
    df_cmb = pd.merge(df, df2, left_index=True, right_index=True)
    print(df_cmb.shape)

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
        if df_corr.loc[idx, col] > 0.9:
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
    ax.set_title('Correlations of genebody mCH \n between known cell types and $de$ $novo$ found clusters')

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

    # tSNE plot as confirmation
    df_celltype = pd.merge(df_cluster, celltype_res, left_on='cluster_ID', right_on='cluster_ID')
    df_plot = pd.merge(df_tsne, df_celltype, left_on='sample', right_on='sample')

    plot_tsne_labels(df_plot, tc='cluster_annotation', 
        title='tSNE of human MB_v1, MB_EA, and MB_EB samples \n (cluster cell type annotation)', 
        figsize=(8,8), legend_mode=1,
        output=output_tsne, show=False
        )

if __name__ == '__main__':

    log = create_logger()
    files = (glob.glob(os.path.join(OUTPUT_DIR, 'test_*p.tsv')) + 
        glob.glob(os.path.join(OUTPUT_DIR, 'test_kmeans_*k.tsv')))
    for i, file in enumerate(files):
        log.info('Processing file {}/{}'.format(i+1, len(files)))
        cluster_fname = os.path.basename(file)
       	input_f = './data/gene_level/genebody_mCH_human_combined_summary.tsv'
        auto_annotate(cluster_fname, input_f, OUTPUT_DIR)
        plt.close('all')

#!/usr/bin/env python3
"""
louvain clustering from kNN graph

2 model parameters: # of PCs and # of nearest neighbors

"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import louvain
import igraph as ig

from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors

import time
import argparse
import logging

from snmcseq_utils import create_logger
from snmcseq_utils import plot_tsne_labels


# jaccard index 
def compute_jaccard(x, y):
    """
    Compute the jaccard index for 2 boolean arrays
    """
    x = np.asarray(x).astype(np.bool)
    y = np.asarray(y).astype(np.bool)
    assert x.shape == y.shape

    intersection = np.logical_and(x, y)
    union = np.logical_or(x, y)
    
    return intersection.sum()/float(union.sum()) 
    
def compute_jaccard_weights_old(X):
    """
    Given an (unweighted) ajacency matrix, assign jaccard index as weights
    
    X is unsymmetric and 0-1 valued 
    """
    X = np.asarray(X)
    Y = np.zeros(X.shape)
    
    # non-zero element index (i, j)
    for i, j in zip(*np.nonzero(X)): 
        Y[i,j] = compute_jaccard(X[i, :], X[j, :])
    return Y

def compute_jaccard_weights(X, option='DIRECTED'):
    """
    Given an (unweighted) ajacency matrix, assign jaccard index as weights
    
    X has to be 0-1 valued 
    """
    X = np.asarray(X)
    ni, nj = X.shape
    assert ni == nj
    
    Y = np.dot(X, X.T)/(ni - np.dot((X-1), (X-1).T))
    
    if option == 'DIRECTED':
        # only set weight in places where X are nonzero
        Y = Y*X
    elif option == 'UNDIRECTED':
        Y = Y*((X+X.T)>0)
    else:
        raise ValueError('choose DIRECTED or UNDIRECTED')
    return Y

def gen_knn_graph(X, k, option='DIRECTED'):
    """
    Generate knn graph from data matrix (n_obs, n_features)
    Edges are weighted by jaccard index between 2 nodes

    X: array like (n_obs, n_features) 
    k: int, number of nearest neighbors
    option: directed or indirected 
    """

    # get 0-1 adjacency matrix of kNN 
    knn = NearestNeighbors(n_neighbors=k, metric='euclidean').fit(X)
    g_knn = knn.kneighbors_graph(X, mode='connectivity')
    g_knn = g_knn.toarray()

    # get weights 
    gw_knn = compute_jaccard_weights(g_knn, option=option)

    return gw_knn

# louvain clustering
def louvain_clustering(adj_mtx, index, option='DIRECTED'):
    """
    a wrap-up function for louvain clustering given a (weighted) adjacency matrix
    """
    G = ig.Graph.Weighted_Adjacency(adj_mtx.tolist(), mode=option)
    louvain.set_rng_seed(1)
    partition1 = louvain.find_partition(G, louvain.ModularityVertexPartition, weights=G.es["weight"])

    labels = [0]*index.shape[0] 
    for i, cluster in enumerate(partition1):
        for element in cluster:
            labels[element] = 'cluster_' + str(i+1)
        
    df_res = pd.DataFrame(index=[idx[:-len('_mcc')] for idx in index])
    df_res['cluster_ID'] = labels 
    df_res = df_res.rename_axis('sample', inplace=True)
    return df_res

# argparser
def create_parser():
    """
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", 
        help="input file (mcc file)", required=True)
    parser.add_argument("-o", "--output", 
        help="output file (clustering results)", required=True)

    parser.add_argument("-n", "--numpc", 
        help="number of PCs", type=int, default=50)
    parser.add_argument("-k", "--knn", 
        help="number of nearest neighbors", type=int, default=30)
    parser.add_argument("-op", "--option", 
        help="DIRECTED or UNDIRECTED", default='UNDIRECTED')
    return parser

if __name__ == '__main__':

    ti = time.time()
    log = create_logger()
    parser = create_parser()
    args = parser.parse_args()

    data_f = args.input
    output_f = args.output
    n = args.numpc
    k = args.knn
    option = args.option

    # import real data
    log.info('Importing data...')
    # mCC
    df = pd.read_table(data_f, index_col=['chr', 'bin'], dtype={'chr': object})
    df = df.T
    print('Input shape (n_obs, n_features): %s' % (df.shape, ))

    # pca
    log.info('PCA...')
    pca = PCA(n_components=n).fit(df.values)
    pcs = pca.transform(df.values)

    # build a jaccard-index weighted k nearest neighbor graph
    log.info('Generating kNN grpah...')
    adj_mtx = gen_knn_graph(pcs, k, option='DIRECTED')

    # run louvain jaccard clustering
    log.info('Louvain clustering...')
    df_res = louvain_clustering(adj_mtx, df.index, option='DIRECTED')

    # output
    log.info('Outputing clustering results to: %s' % output_f)
    df_res.to_csv(output_f,
                 sep='\t', na_rep='NA', header=True, index=True)

    tf = time.time()
    log.info("Running time: %.2f second s" % (tf-ti))


    # # plot result
    # import tSNE coords
    # tsne_f = '/cndd/fangming/snmcseq_dev/data/tsne/tsne_perp30_binc_mCH_human_combined_100000_summary_nmcc_v3.tsv'
    # # tsne_f = '/cndd/fangming/snmcseq_dev/data/tsne/tsne_perp30_binc_mCG_human_combined_100000_summary_nmcc_v3.tsv'
    # # tsne_f = '/cndd/fangming/snmcseq_dev/data/tsne/tsne_perp30_binc_mCHmCG_human_combined_100000_summary_nmcc_v3.tsv'
    # df_tsne = pd.read_table(tsne_f, index_col='sample')
    # print('Input tSNE file shape: %s' % (df_tsne.shape, ))

    # df_plot = pd.merge(df_tsne, df_res, left_index=True, right_index=True)
    # plot_tsne_labels(df_plot, tc='cluster_ID', legend_mode=1)

    # plot_tsne_labels(df_plot, tc='cluster_ID', 
    #     title='tSNE of human MB_v1, MB_EA, and MB_EB samples \n (clusters generated by Louvain method on kNN graph)', 
    #     figsize=(8,8), legend_mode=1,
    #     output='/cndd/fangming/snmcseq_dev/results/tsne_clusters/tsne_cluster_MB_v1_MB_EA_MB_EB_cluster_v1.pdf'
    #     )



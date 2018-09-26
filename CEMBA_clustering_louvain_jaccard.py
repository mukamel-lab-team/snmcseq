#!/usr/bin/env python3
"""
louvain clustering from kNN graph

2 model parameters: # of PCs and # of nearest neighbors

"""

from __init__ import *

from collections import OrderedDict
import argparse

from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
import louvain
import igraph as ig
from scipy import sparse

from snmcseq_utils import create_logger
from snmcseq_utils import plot_tsne_labels


def compute_jaccard_weights_v2(X, k):
    """compute jaccard index on a knn graph
    Arguments: 
        X (unweighted) kNN ajacency matrix (each row Xi* gives the kNNs of cell i) 
        X has to be 0-1 valued 
        k (number of nearest neighbors) 
        
    output: numpy matrix Y
    """
    X = sparse.csr_matrix(X)
    ni, nj = X.shape
    assert ni == nj
    
    tmp = X.dot(X.T)
    Y = X.multiply(tmp/(2*k - tmp.todense()))    
    return Y

def compute_jaccard_weights(X, option='DIRECTED'):
    """
    Given an (unweighted) ajacency matrix, assign jaccard index as weights
    
    X has to be 0-1 valued 
    """
    X = np.asarray(X).astype(float)
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
    # get weights 

    knn = NearestNeighbors(n_neighbors=k, metric='euclidean').fit(X)
    g_knn = knn.kneighbors_graph(X, mode='connectivity')
    g_knn = g_knn.toarray()
    gw_knn = compute_jaccard_weights(g_knn, option=option)
    
    # updated 05/15/2018 Fangming
    # gw_knn = compute_jaccard_weights_v2(g_knn, k=k).todense()

    return gw_knn

# louvain clustering
def louvain_clustering(adj_mtx, index, option='DIRECTED', sample_column_suffix=None):
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

    if sample_column_suffix: 
        df_res = pd.DataFrame(index=[idx[:-len(sample_column_suffix)] for idx in index])
    else:
        df_res = pd.DataFrame(index=index)
    df_res['cluster_ID'] = labels 
    df_res = df_res.rename_axis('sample', inplace=True)
    return df_res

def louvain_jaccard(df, n_pc=50, k=30, sub_ncells=None, output_file=None, sample_column_suffix='_mcc', whiten=False):
    """louvain jaccard clustering from feature matrix
    df is a gene-by-cell or bin-by-cell matrix 
    cells column has the pattern of "_mcc$"

    """
    ti = time.time()
    df = df.T

    # subsample cells
    if sub_ncells:
        df = df.sample(n=sub_ncells, random_state=1)

    logging.info('Begin louvain jaccard clustering\nInput shape (n_obs, n_features): {}'.format(df.shape))

    # pca
    pca = PCA(n_components=n_pc, whiten=whiten).fit(df.values)
    pcs = pca.transform(df.values)

    # build a jaccard-index weighted k nearest neighbor graph
    adj_mtx = gen_knn_graph(pcs, k, option='DIRECTED')

    # run louvain jaccard clustering
    tii = time.time()
    df_res = louvain_clustering(adj_mtx, df.index, option='DIRECTED', sample_column_suffix=sample_column_suffix)

    # number of clusters
    nclst = np.unique(df_res.cluster_ID.values).shape[0]

    # total time
    tf = time.time()
    
    # summary 
    summary = OrderedDict({'time': tf-ti, 'time_clst': tf-tii, 'nclst': nclst, 'n_cells': df.shape[0],
                           'n_pc': n_pc, 'k': k})
    logging.info('clustering summary: {}'.format(summary))

    # output
    if output_file:
        logging.info('Outputing clustering results to: {}'.format(output_file))
        df_res.to_csv(output_file,
                     sep='\t', na_rep='NA', header=True, index=True)
        # summary
        keys = list(summary.keys()) 
        values = [str(value) for value in list(summary.values())] 
        with open(output_file+'.summary', 'w') as file:
            file.write('\t'.join(keys)+'\n')
            file.write('\t'.join(values)+'\n') 

    return df_res, summary

def plot_louvain_jaccard(df_cluster, df_tsne):
    """
    """
    pass

def run_louvain_CEMBA(ens, ks=K_NN, n_pc=N_PC):
    """
    """
    ens_path = os.path.join(PATH_ENSEMBLES, ens)
    nmcc_files = sorted(glob.glob(os.path.join(ens_path, 'binc/binc_*_nmcc_{}.tsv'.format(ens)))) 

    if not os.path.isdir(os.path.join(ens_path, 'cluster')):
        os.makedirs(os.path.join(ens_path, 'cluster'))
    if not os.path.isdir(os.path.join(ens_path, 'plots')):
        os.makedirs(os.path.join(ens_path, 'plots'))

    for nmcc_file in nmcc_files:
        nmcc_basename = os.path.basename(nmcc_file) 
        df = pd.read_table(nmcc_file, index_col=['chr', 'bin'], dtype={'chr': object, 'bin': int})
        for k in ks:
            output_file = os.path.join(ens_path, 'cluster/cluster_lv_npc{}_k{}_{}'.format(n_pc, k, nmcc_basename))
            df_cluster = louvain_jaccard(df, n_pc=n_pc, k=k, sub_ncells=None, output_file=output_file)
    return


if __name__ == '__main__':

    ti = time.time()
    log = create_logger()
    log.info("Begin louvain clustering...")

    enss = ['Ens1', 'Ens2', 'Ens3', 'Ens4']
    enss = ['Ens0']
    for ens in enss:    
        run_louvain_CEMBA(ens)

    tf = time.time()
    log.info("Total louvain clustering running time: {} seconds".format(tf-ti))



#!/usr/bin/env python3
"""
louvain clustering from kNN graph
2 model parameters: # of PCs and # of nearest neighbors

updated from CEMBA_clustering_louvain_jaccard.py for large number of cells
"""

from __init__ import *
# from collections import OrderedDict
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
import louvain
import igraph as ig
from scipy import sparse

from snmcseq_utils import create_logger

def gen_knn_annoy_train_test(X_train, X_test, k, form='list', metric='euclidean', n_trees=10, search_k=-1, verbose=True):
    """X is expected to have low feature dimensions (n_obs, n_features) with (n_features <= 50)
    For each row in X_test, find k nearest neighbors in X_train
    """
    from annoy import AnnoyIndex
    
    ti = time.time()
    
    n_obs, n_f = X_train.shape
    n_obs_test, n_f_test = X_test.shape
    assert n_f == n_f_test 
    
    t = AnnoyIndex(n_f, metric=metric)  # Length of item vector that will be indexed
    for i, X_row in enumerate(X_train):
        t.add_item(i, X_row)

    t.build(n_trees) # 10 trees
    knn = [0]*(n_obs_test)
    for i, X_row_test in enumerate(X_test):
        knn[i] = t.get_nns_by_vector(X_row_test, k, search_k=search_k) # will find the 1000 nearest neighbors
    knn = np.array(knn)
   
    if verbose:
        print("Time used to get kNN {}".format(time.time()-ti))
    
    if form == 'adj':
        # row col 1 
        row_inds = np.repeat(np.arange(n_obs_test), k)
        col_inds = np.ravel(knn)
        data = [1]*len(row_inds)
        return sparse.coo_matrix((data, (row_inds, col_inds)), shape=(n_obs_test, n_obs))
    elif form == 'list':  #
        return knn
    else:
        raise ValueError("Choose from 'adj' and 'list'")

def gen_knn_annoy(X, k, form='list', metric='euclidean', n_trees=10, search_k=-1, verbose=True):
    """X is expected to have low feature dimensions (n_obs, n_features) with (n_features <= 50)
    """
    from annoy import AnnoyIndex

    ti = time.time()

    n_obs, n_f = X.shape
    k = min(k, n_obs)
    if k == n_obs:
        print("Actual k: {}".format(n_obs))
     
    t = AnnoyIndex(n_f, metric=metric)  # Length of item vector that will be indexed
    for i, X_row in enumerate(X):
        t.add_item(i, X_row)

    t.build(n_trees) # 10 trees
    knn = [0]*(n_obs)
    for i in range(n_obs):
        knn[i] = t.get_nns_by_item(i, k, search_k=search_k) # will find the 1000 nearest neighbors
    knn = np.array(knn)
   
    if verbose:
        print("Time used to get kNN {}".format(time.time()-ti))
    
    if form == 'adj':
        # row col 1 
        row_inds = np.repeat(np.arange(n_obs), k)
        col_inds = np.ravel(knn)
        data = [1]*len(row_inds)
        return sparse.coo_matrix((data, (row_inds, col_inds)), shape=(n_obs, n_obs))
    elif form == 'list':  #
        return knn
    else:
        raise ValueError("Choose from 'adj' and 'list'")

def gen_knn(pcX, k, form='adj', metric='euclidean', verbose=True): 
    """Generate kNN matrix from a pcX (n_obs, n_feature) matrix
    
    """
    ti = time.time()

    n_obs, n_f = X.shape
    k = min(k, n_obs)
    
    knn = NearestNeighbors(n_neighbors=k, metric=metric).fit(pcX)
    
    if form == 'adj':
        g_knn = knn.kneighbors_graph(pcX, mode='connectivity')
        if verbose:
            print("Time spent on generate kNN graph: {}".format(time.time()-ti))
        return g_knn
            
    elif form == 'list':
        dists, inds = knn.kneighbors(pcX)
        if verbose:
            print("Time spent on generate kNN graph: {}".format(time.time()-ti))
        return (dists, inds) 
    
def compute_jaccard_weights_from_knn(X):
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
    
    k = X[0, :].sum() # number of neighbors
    
    Y = X.dot(X.T)
    # Y = X.multiply(tmp/(2*k - tmp.todense()))    
    Y.data = Y.data/(2*k - Y.data)
    
    return Y 

def adjacency_to_igraph(adj_mtx, weighted=False):
    """
    Converts an adjacency matrix to an igraph object
    
    Args:
        adj_mtx (sparse matrix): Adjacency matrix
        directed (bool): If graph should be directed
    
    Returns:
        G (igraph object): igraph object of adjacency matrix
    
    Uses code from:
        https://github.com/igraph/python-igraph/issues/168
        https://stackoverflow.com/questions/29655111

    Author:
        Wayne Doyle 
        (Fangming Xie modified) 
    """
    nrow, ncol = adj_mtx.shape
    if nrow != ncol:
        raise ValueError('Adjacency matrix should be a square matrix')
    vcount = nrow
    sources, targets = adj_mtx.nonzero()
    edgelist = list(zip(sources.tolist(), targets.tolist()))
    G = ig.Graph(n=vcount, edges=edgelist, directed=True)
    if weighted:
        G.es['weight'] = adj_mtx.data
    return G

def louvain_lite(G, cell_list, weighted=False, verbose=True):
    """
    weighted=False is 10x faster than True
    """
    ti = time.time()
        
    if weighted:
        partition1 = louvain.find_partition(G, louvain.ModularityVertexPartition, 
                                            weights=G.es["weight"]
                                           )
    else:
        partition1 = louvain.find_partition(G, louvain.ModularityVertexPartition) 
        
    labels = [0]*(len(cell_list)) 
    for i, cluster in enumerate(partition1):
        for element in cluster:
            labels[element] = i+1

    df_res = pd.DataFrame(index=cell_list)
    df_res['cluster'] = labels 
    df_res = df_res.rename_axis('sample', inplace=False)
    
    if verbose:
        print("Time spent on louvain clustering: {}".format(time.time()-ti))
        
    return df_res

def clustering_routine(X, cell_list, k, metric='euclidean', option='plain', n_trees=10, search_k=-1):
    """
    X is a (n_obs, n_feature) matrix, n_feature <=50 is recommended
    option: {'plain', 'jaccard', ...}
    """
    if option == 'plain':
        g_knn = gen_knn_annoy(X, k, form='adj', metric=metric, 
                              n_trees=n_trees, search_k=search_k, verbose=True)
        G = adjacency_to_igraph(g_knn, weighted=False)
        df_res = louvain_lite(G, cell_list, weighted=False)
        
    elif option == 'jaccard':
        g_knn = gen_knn_annoy(X, k, form='adj', metric=metric, 
                              n_trees=n_trees, search_k=search_k, verbose=True)
        gw_knn = compute_jaccard_weights_from_knn(g_knn)
        G = adjacency_to_igraph(gw_knn, weighted=True)
        df_res = louvain_lite(G, cell_list, weighted=True)
    else:
        raise ValueError('Choose from "plain" and "jaccard"')
    
    return df_res

def clustering_routine_old(X, cell_list, k, metric='euclidean', option='plain'):
    """
    X is a (n_obs, n_feature) matrix, n_feature <=50 is recommended
    option: {'plain', 'jaccard', ...}
    """
    if option == 'plain':
        g_knn = gen_knn(X, k, form='adj', metric=metric, verbose=True)
        G = adjacency_to_igraph(g_knn, weighted=False)
        df_res = louvain_lite(G, cell_list, weighted=False)
        
    elif option == 'jaccard':
        g_knn = gen_knn(X, k, form='adj', metric=metric, verbose=True)
        gw_knn = compute_jaccard_weights_from_knn(g_knn)
        G = adjacency_to_igraph(gw_knn, weighted=True)
        df_res = louvain_lite(G, cell_list, weighted=True)
    else:
        raise ValueError('Choose from "plain" and "jaccard"')
    
    return df_res




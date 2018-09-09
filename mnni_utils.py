#!/usr/bin/env python3

from __init__ import *
from scipy import sparse
from sklearn.decomposition import PCA


def smooth_in_modality(counts_matrix, norm_counts_matrix, k, ka, npc=100, sigma=1.0, p=0.9, r=None):
    """Smooth a data matrix
    
    Arguments:
        - counts_matrix (pandas dataframe, feature by cell)
        - norm_counts_matrix (pandas dataframe, feature by cell) log10(CPM+1)
        - k (number of nearest neighbors)
    Return:
        - smoothed cells_matrix (pandas dataframe)
        - markov affinity matrix
    """
    from sklearn.neighbors import NearestNeighbors
    
    assert counts_matrix.shape[1] == norm_counts_matrix.shape[1] 

    c = norm_counts_matrix.columns.values
    N = len(c)

    # reduce dimension
    pca = PCA(n_components=npc, random_state=1) 
    pcs = pca.fit_transform(norm_counts_matrix.T)

    # get k nearest neighbor distances 
    knn = NearestNeighbors(n_neighbors=k, metric='euclidean').fit(pcs)
    dists, inds = knn.kneighbors(pcs)

    # remove itself
    dists = dists[:, 1:]
    inds = inds[:, 1:]

    # normalize by ka's distance 
    dists = (dists/(dists[:, ka].reshape(-1, 1)))

    # gaussian kernel
    adjs = np.exp(-((dists**2)/(sigma**2))) 

    # construct a sparse matrix 
    cols = np.ravel(inds)
    rows = np.repeat(np.arange(N), k-1) # remove itself
    vals = np.ravel(adjs)
    A = sparse.csr_matrix((vals, (rows, cols)), shape=(N, N))

    # Symmetrize A
    A = A + A.T

    # normalization (A is now a weight matrix excluding itself)
    A = sparse.csr_matrix(A/A.sum(axis=1))

    # include itself
    eye = sparse.identity(N)
    if p:
        A = p*eye + (1-p)*A
    elif r:
        A = (1.0/(1.0+r))*eye + (r/(1.0+r))*A
    
    # smooth  
    counts_matrix_smoothed = pd.DataFrame((A.dot(counts_matrix.T)).T, 
                                         columns=counts_matrix.columns, index=counts_matrix.index)
    
    return counts_matrix_smoothed, A

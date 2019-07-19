from __init__ import *
from __init__jupyterlab import *
import snmcseq_utils

def select_hvg(gbc_cpm, percentile=30, n_qcut=20, ylim=[]):
    # further select highly variable genes
    # variance/mean
    from sklearn.utils.sparsefuncs import mean_variance_axis
    mean_cpm, var_cpm = mean_variance_axis(gbc_cpm.data.tocsr(), axis=1) 
    vmr_cpm = (var_cpm+1)/(mean_cpm+1)
    
    # select top 30 percentile vmr from each first 9 deciles of CPM
    _x = pd.qcut(pd.Series(mean_cpm), n_qcut, labels=False).to_frame('decile')
    hvgs = []
    fig, ax = plt.subplots()
    for decile, _x_sub in _x.groupby('decile'):
        gene_group = _x_sub.index.values
        mean_cpm_gg = mean_cpm[gene_group]
        vmr_cpm_gg = vmr_cpm[gene_group]
        # genes with top 30% of vmr
        hvg_group = gene_group[vmr_cpm_gg > np.percentile(vmr_cpm_gg, 100-percentile)]

        if decile != n_qcut-1:
            hvgs.append(hvg_group)
            ax.scatter(mean_cpm_gg, vmr_cpm_gg, s=1, c='black', alpha=0.5)
            ax.scatter(mean_cpm[hvg_group], vmr_cpm[hvg_group], s=2)

    hvgs = np.hstack(hvgs)
    
    if ylim:
        ax.set_ylim(ylim)
    plt.show()
    
    return hvgs

def select_hvg_methylation(df_nmcc, percentile=30, n_qcut=20):
    # further select highly variable genes
    # standard deviation 
    
    stds_nmcc = df_nmcc.std(axis=1)
    mean_nmcc = df_nmcc.mean(axis=1)

    # select top 30 percentile vmr from each first 9 deciles of NMCC 
    _x = pd.qcut(mean_nmcc, n_qcut, labels=False).to_frame('decile')
    hvgs = []
    fig, ax = plt.subplots()
    for decile, _x_sub in _x.groupby('decile'):
        gene_group = _x_sub.index.values

        mean_nmcc_gg = mean_nmcc.loc[gene_group]
        stds_nmcc_gg = stds_nmcc.loc[gene_group]
#         print(gene_group.shape, stds_nmcc_gg.shape)
        # genes with top 30% of stds 
        hvg_group = gene_group[stds_nmcc_gg > np.percentile(stds_nmcc_gg, 100-percentile)]

    #     if decile != n_qcut-1:
        hvgs.append(hvg_group)
        ax.scatter(mean_nmcc_gg, stds_nmcc_gg, s=1, c='black', alpha=0.5)
        ax.scatter(mean_nmcc.loc[hvg_group], stds_nmcc.loc[hvg_group], s=2)

    hvgs = np.hstack(hvgs)

    plt.show()
    return hvgs

def filter_genes(gxc_raw, sufficient_cell_coverage=0.01):
    """
    """
    n_gene, n_cell = gxc_raw.data.shape
    gene_cov = (gxc_raw.data > 0).sum(axis=1)
    gene_cov = np.array(gene_cov).squeeze()/n_cell # fraction of cells covered
    cond = gene_cov>sufficient_cell_coverage
    gxc_raw_filtered = GC_matrix(np.array(gxc_raw.gene)[cond],
                                gxc_raw.cell,
                                gxc_raw.data.tocsr()[cond, :],
                               )
    return gxc_raw_filtered

def preproc_rna_cpm_based(gxc_raw, sufficient_cell_coverage=0.01, 
                          hv_percentile=30, hv_ncut=20):
    # select genes expressed in > 1% of cells
    print("Removing low coverage genes...")
    lib_size = np.ravel(gxc_raw.data.sum(axis=0))
    _gxc_tmp = filter_genes(gxc_raw, sufficient_cell_coverage=sufficient_cell_coverage)
    
    # CPM matrix
    print("Getting CPM..")
    gxc_ftr = snmcseq_utils.sparse_logcpm(_gxc_tmp, mode='cpm', lib_size=lib_size)
    del _gxc_tmp

    # select highy variable genes
    print("Getting highly variable genes and logCPM...")
    hvgs = select_hvg(gxc_ftr, percentile=hv_percentile, n_qcut=hv_ncut)
    
    gxc_hvftr = GC_matrix(
                          gxc_ftr.gene[hvgs],
                          gxc_ftr.cell,
                          gxc_ftr.data.tocsr()[hvgs, :],
                          )
    del gxc_ftr
    gxc_hvftr.data.data = np.log10(1+gxc_hvftr.data.data) # very important 
    print("Number of genes: {}".format(len(hvgs)))
    return gxc_hvftr

def preproc_rna_tpm_based(gxc_raw, gene_lengths, impute_length=True, 
                          impute_gene_lengths=True, 
                          sufficient_cell_coverage=0.01, 
                          hv_percentile=30, hv_ncut=20):
    """Gene lengths is a gene length pandas series indexed by gene names
    """
    assert np.all(gxc_raw.gene == gene_lengths.index.values) 
    if impute_gene_lengths:
        print("Imputing gene lengths...")
        gene_lengths = gene_lengths.fillna(gene_lengths.mean())
    lib_size = np.ravel(gxc_raw.data.sum(axis=0))

    # select genes expressed in > 1% of cells
    print("Removing low coverage genes...")
    _gxc_tmp = filter_genes(gxc_raw, sufficient_cell_coverage=sufficient_cell_coverage)
    
    # CPM matrix
    print("Getting CPM..")
    gxc_ftr = snmcseq_utils.sparse_logcpm(_gxc_tmp, mode='cpm', lib_size=lib_size)
    del _gxc_tmp

    # select highy variable genes
    print("Getting highly variable genes and logCPM...")
    hvgs = select_hvg(gxc_ftr, percentile=hv_percentile, n_qcut=hv_ncut) # index in gxc_ftr
    hvgs_genes = gxc_ftr.gene[hvgs]
    del gxc_ftr

    # TPM matrix from gxc_raw
    print("Getting logTPM...")
    gxc_logtpm = snmcseq_utils.sparse_logtpm(gxc_raw, gene_lengths)
    hvgs_idx = snmcseq_utils.get_index_from_array(gxc_logtpm.gene, hvgs_genes)

    # Trim logTPM matrix
    print("Trim logTPM matrix...")
    gxc_hvftr = GC_matrix(
                          gxc_logtpm.gene[hvgs],
                          gxc_logtpm.cell,
                          gxc_logtpm.data.tocsr()[hvgs, :],
                          )
    print("Number of genes: {}".format(len(hvgs)))
    return gxc_hvftr

def preproc_methylation():
    pass
#!/usr/bin/env python3

from __init__ import *

from scipy.stats.mstats import zscore
from scipy.stats import entropy 

import snmcseq_utils
from CEMBA_update_mysql import connect_sql

def counts_to_nmcc(cells_mc_c, 
                   df_cells, clst_col, global_mcc_col,
                   base_call_cutoff=20, sufficient_coverage_fraction=0.9, 
                   normalize=True, 
                   suffix=False, 
                   skip_cell_level=False,
                  ):
    """Given count matrix and global mcc levels, prepare cells_nmcc and clusters_nmcc
    
    Args:
        - cells_mc_c: dataframe, indexed by feature, columns are cell_name 
        - df_cells: dataframe, indexed by cell_name, columns are ['cluster', 'global_mcc']
            cluster assignment, global_mcc 
        
        - base_call_cutoff
        - sufficient_coverage_fraction,
        - suffix
        
        - skip_cell_level: default False, if True, it skips cells_mcc=None
    
    Return:
        - cells_nmcc: dataframe, indexed by features (selected), columns are cell_name
        - clusters_nmcc: dataframe, indexed by features (selected), columns are cell_name
    """
    
    # check cell name matched
    assert len(df_cells.index) == len(cells_mc_c.columns)/2
    
    # get mcc (apply loose threshold)
    if not skip_cell_level:
        cells_mcc = snmcseq_utils.get_mcc(cells_mc_c, 
                            base_call_cutoff=base_call_cutoff, 
                            sufficient_coverage_fraction=sufficient_coverage_fraction, 
                            suffix=suffix)
    else:
        cells_mcc = None
    
    # cluster_mcc
    clusters_mc_c = snmcseq_utils.get_clusters_mc_c_worker(df_cells, cells_mc_c, clst_col)
    clusters_mcc = snmcseq_utils.get_mcc(clusters_mc_c, suffix=suffix)
    
    # str to int if cluster names are pure number
    try:
        cols = [int(col) for col in clusters_mcc.columns]
    except:
        cols = clusters_mcc.columns 
    clusters_mcc.columns = cols 

    if normalize:
        # get nmcc
        if not skip_cell_level:
            cells_nmcc = cells_mcc.divide(df_cells[global_mcc_col], axis=1)
        else:
            cells_nmcc = None
        
        clusters_global_mcc = df_cells.groupby(clst_col).median()[global_mcc_col]
        clusters_nmcc = clusters_mcc.divide(clusters_global_mcc, axis=1)
        
        return cells_nmcc, clusters_nmcc
    else:
        return cells_mcc, clusters_mcc


def find_marker_genes(cells_nmcc, clusters_nmcc, df_cells, df_genes, 
                        ntop=30, p_std=0.75, p = 0.5, p_putative=0.05, bins=np.arange(0, 3, 0.05)):
    """Given cluster assignment, gene annotation and data matrix, 
    compute top marker genes for each cluster.
    Metrics: KL divergence and difference in mean/std
    
    This could in principle be extended to compute any "marker features", but coverage might be an issue 
    for short genomic region as it requires the feature to have good coverage at cell level.
    
    Args:
        - cells_nmcc matrix: 
            A dataframe of nmcc values of gene by cell matrix 
            index: gene_id; columns: [cell_name1, cell_name2, ...]
        
        - clusters_nmcc matrix: 
            A dataframe of nmcc values of gene by cluster matrix 
            index: gene_id; columns: [cluster1, cluster2, ...]
        
        - df_cells: cluster assignment of cells  
            A dataframe with cell name as index and cluster assignment as a column.
            index: cell_name; columns: ['cluster']
            
        - df_genes: gene_annotation 
            A dataframe with gene_id as index and gene_name as a column.
            index: gene_id; columns: ['gene_name']
        
        - ntop: default 30
            number of top marker genes
        - p_std: default 0.75
            threshold percentile of standard deviation above which mean difference is penalized 
        - p: default 0.5
            fraction of mean_dev (overall score = p*mean_dev + (1-p)*kl_divergence)
        - p_putative: default 0.05
            fraction of genes go into further consideration
        - bins: default = np.arange(0, 3, 0.05)
            bins used for kl divergence calculation
            
    Return:
        - A dataframe of top marker genes for each cluster/column
    """

    # check
    # gene_id as index, 'gene_name'
    assert 'gene_name' in df_genes.columns
    # 'cluster'
    assert 'cluster' in df_cells.columns
    # cell names match
    assert sorted(cells_nmcc.columns.tolist()) == sorted(df_cells.index.tolist())
    # cluster names match
    assert sorted(clusters_nmcc.columns.tolist()) == sorted(df_cells['cluster'].unique().tolist())
    
    markers_all = pd.DataFrame() 
    
    # zscore by cluster, quantile by genes
    data_pct = clusters_nmcc.apply(lambda x: zscore(x, ddof=1), axis=1).rank(pct=True, axis=0)

    for col in clusters_nmcc.columns:
        print(col, end=',')

        # putative markers (initial filtering)
        marker_ids = (1 - data_pct[(data_pct[col] < p_putative)]).index
        
        # dropna is necessary to exclude genes that don't have coverage at cell level 
        data_markers = cells_nmcc.loc[marker_ids, :].dropna()
        data_clst_markers = clusters_nmcc.loc[marker_ids, :].dropna()

        marker_ids = np.intersect1d(data_markers.index, data_clst_markers.index)
        markers = df_genes.loc[marker_ids, :].copy() 
        data_markers = cells_nmcc.loc[marker_ids, :]
        data_clst_markers = clusters_nmcc.loc[marker_ids, :]

        # mean dev 
        stds = data_clst_markers.std(axis=1)
        bar_std = np.percentile(stds, 100*p_std)
        stds_masked = np.maximum(stds, bar_std)

        mean_dev = (data_clst_markers[col] - data_clst_markers.mean(axis=1)).divide(stds_masked, axis=0)

        # markers ranking by KL divergence 
        # px
        pxs = data_markers.apply(
            lambda x: (np.histogram(x, bins=bins)), axis=1).apply(
            lambda x: (x[0]+1)/((x[0]+1).sum()))

        # py
        # cells NOT in the cluster
        cells = df_cells[df_cells['cluster']!=col].index.values
        pys = data_markers.loc[:, cells].apply(
            lambda x: (np.histogram(x, bins=bins)), axis=1).apply(
            lambda x: (x[0]+1)/((x[0]+1).sum()))

        # KL
        ps = pd.concat((pxs, pys), axis=1)
        kls = ps.apply(lambda x: entropy(x[0], x[1]),  axis=1)


        ### markers rank by overall zscore
        markers['zscore_mean_dev'] = -zscore(mean_dev.values) # aim high 
        markers['zscore_kl'] = zscore(kls.values) # aim high

        markers['overall_score'] = (p*markers['zscore_mean_dev'] 
                                    + (1-p)*markers['zscore_kl']
                                   )
        markers = markers.sort_values('overall_score', ascending=False) # aim high

        # fix bug: head(ntop) might be less than ntop
        res = markers.head(ntop)['gene_name'].values
        if ntop - len(res) > 0:
	        res = np.hstack((res, np.repeat(np.nan, (ntop-len(res)))))
        markers_all[col] = res  

    print('\nDone.')
    return markers_all




def find_marker_genes_CEMBA(ens, context='CH', clsts='auto', 
                            base_call_cutoff=20, sufficient_coverage_fraction=0.90,
                            ntop=30, p_std=0.75, p=0.5, p_putative=0.05, bins=np.arange(0, 3, 0.05),
                           ):
    """A wrap-up function of `find_marker_genes` for CEMBA project
    
    It reads info as needed from CEMBA mySQL and CEMBA directory, 
    and outputs marker gene tables to ens/cluster directory
    in the format of $clst_filename.$context_marker_genes.
    
    Args:
        - ens: ensemble_id, e.g. "Ens1"
        - context: default 'CH'
        - clsts" default 'auto' (get all clusters from CEMBA_mySQL)
        
        - ntop: number of top marker genes
        - other parameters -> function find_marker_genes for more details 
    Return:
        None
    """
    logging.info("Marker genes: {}".format(ens))

    engine = connect_sql(DATABASE)
    
    # get cells info
    ens_path = os.path.join(PATH_ENSEMBLES, ens)
    sql = """SELECT cell_name, dataset, global_m{1}, {0}.* FROM cells
            JOIN {0} ON cells.cell_id = {0}.cell_id""".format(ens, context)
    df_cells = pd.read_sql(sql, engine, index_col='cell_name')
    
    # get gene annotation
    sql = '''SELECT * FROM genes'''
    df_genes = pd.read_sql(sql, engine, index_col='gene_id')
    
    # get cells_mc_c (might be time consuming)
    cells_mc_c = snmcseq_utils.pull_genebody_mc_c(ens, context=context)
    logging.info("Got cells_mc_c matrix: {}".format(cells_mc_c.shape))
    
    # get cells_nmcc, clusters_nmcc
    # all clusterings available in the ens
    cluster_cols = df_cells.filter(regex='^cluster').columns.values
    for i, cluster_col in enumerate(cluster_cols):
        logging.info('Looking for marker genes: {} ({}/{})'.format(cluster_col, i+1, len(cluster_cols)))
        if i == 0:
            skip_cell_level = False
            # get mcc (apply loose threshold)
            cells_nmcc, clusters_nmcc = counts_to_nmcc(cells_mc_c, 
                           df_cells, clst_col=cluster_col, global_mcc_col='global_m{}'.format(context),
                           base_call_cutoff=base_call_cutoff, 
                        sufficient_coverage_fraction=sufficient_coverage_fraction, 
                           normalize=True, 
                           suffix=False, 
                        skip_cell_level=skip_cell_level, 
                        )
        else:
            skip_cell_level = True
            # get mcc (apply loose threshold)
            dis, clusters_nmcc = counts_to_nmcc(cells_mc_c, 
                           df_cells, clst_col=cluster_col, global_mcc_col='global_m{}'.format(context),
                           base_call_cutoff=base_call_cutoff, 
                        sufficient_coverage_fraction=sufficient_coverage_fraction, 
                           normalize=True, 
                           suffix=False, 
                        skip_cell_level=skip_cell_level, 
                        )

        # find marker genes
        res = find_marker_genes(cells_nmcc, clusters_nmcc, df_cells.rename(columns={cluster_col: 'cluster'}), df_genes, 
                            ntop=ntop, p_std=p_std, p=p, p_putative=p_putative, bins=bins)

        # output to file
        dis, clst_context, method, npc, k = cluster_col.split('_')
        output_file = 'cluster_{}_{}_{}_binc_{}_100000_nmcc_{}.tsv.{}_marker_genes'.format(
                        method, npc, k, clst_context, ens, context)
        output_file = os.path.join(PATH_ENSEMBLES, ens, 'cluster', output_file)
        res.to_csv(output_file, sep='\t', header=True, index=False, na_rep='NA')
        logging.info('Saved to {}'.format(output_file))

    return

if __name__ == '__main__':

	enss = ['Ens{}'.format(i) for i in np.arange(11, 51, 1)]
	# enss = ['Ens10', 'Ens51']
	context = 'CH'

	log = snmcseq_utils.create_logger()
	for ens in enss:
		find_marker_genes_CEMBA(ens, context=context, clsts='auto', p_putative=0.10) 

#!/usr/bin/env python3

import sqlalchemy as sa
from scipy import sparse
from collections import namedtuple

from __init__ import *
from snmcseq_utils import create_logger
from snmcseq_utils import cd 
from CEMBA_update_mysql import insert_into
from CEMBA_update_mysql import connect_sql 
from CEMBA_update_mysql import gene_id_to_table_name 
from CEMBA_update_mysql import upload_to_datasets


# from snmcseq_utils import get_mouse_chromosomes
# from snmcseq_utils import compute_global_mC 
# from CEMBA_run_tsne import run_tsne_CEMBA

# def gene_id_to_name(gene_id, df_genes):
#     """df_genes
#     """
#     return df_genes.loc[gene_id, 'gene_name']

# def gene_name_to_id(gene_name, df_genes):
#     """df_genes
#     """
#     return df_genes[df_genes['gene_name'] == gene_name].index.values[0]

GC_matrix = namedtuple('GC_matrix', ['gene', 'cell', 'data'])
def sparse_logcpm(gc_matrix):
    """
    """
    lib_size_inv = sparse.diags(np.ravel(1.0/gc_matrix.data.sum(axis=0)))
    logcpm = (gc_matrix.data).dot(lib_size_inv*1e6).tocoo()
    logcpm.data = np.log10(logcpm.data + 1)

    gc_logcpm = GC_matrix(
        gc_matrix.gene, 
        gc_matrix.cell, 
        logcpm,
    )
    
    return gc_logcpm

def upload_dataset_worker(dataset, gc_matrix_raw, database=DATABASE_ATAC):
    """Given a gene-by-cell count matrix, calculate log10(CPM+1), upload it to CEMBA_ATAC
    """
    engine = connect_sql(database)

    # logcpm
    gc_logcpm = sparse_logcpm(gc_matrix_raw)

    # upload_to_datasets
    upload_to_datasets(dataset, database=database, strict=False)

    # upload to cells
    df_cells = pd.DataFrame()
    df_cells['cell_name'] = gc_logcpm.cell 
    df_cells['dataset'] = dataset 
    insert_into(engine, 'cells', df_cells, ignore=True, verbose='True')

    # upload to genes
    sql = """SELECT * FROM cells"""
    df_cells = pd.read_sql(sql, engine)
    metadata = sa.MetaData(engine)

    for i, gene in enumerate(gc_logcpm.gene):
        # if i > 10:
        #     break
        
        data_gene = gc_logcpm.data.getrow(i).tocoo()
        gene_table_name = gene_id_to_table_name(gene)

        if (i%100 == 0):
            logging.info("Progress on genes: {} {}".format(i+1, gene_table_name))
        
        df_gene = pd.DataFrame()
        df_gene['cell_name'] = [gc_logcpm.cell[col] for col in data_gene.col]
        df_gene['normalized_counts'] = data_gene.data
        
        df_gene = pd.merge(df_gene, df_cells, on='cell_name')
        df_gene = df_gene[['cell_id', 'normalized_counts']]

        # print(gene, df_gene)
        if not df_gene.empty:
            insert_into(engine, gene_table_name, df_gene, ignore=True, verbose=False)

    return

def read_sparse_atac_matrix(fdata, frow, fcol, dataset):
    """
    Example:
        dataset = 'CEMBA_1B_180118'
        fdir = '/cndd/Public_Datasets/CEMBA/snATACSeq/Datasets/1B/CEMBA180118_1B_rep1/counts/genebody/'
        fdata = fdir + 'CEMBA180118_1B_rep1.genebody.count_matrix.npz'
        frow = fdir + 'CEMBA180118_1B_rep1.genebody.count_matrix.row.index' # gene
        fcol = fdir + 'CEMBA180118_1B_rep1.genebody.count_matrix.col.index' # col

    """

    gene_ids = pd.read_table(frow, header=None)[0].values
    cells = pd.read_table(fcol, header=None)[0].values
    cells = cells + '_' + dataset 
    data = sparse.load_npz(fdata)

    gc_matrix_raw = GC_matrix(
        gene_ids, 
        cells, 
        data,
    )
    assert gc_matrix_raw.data.shape[0] == len(gc_matrix_raw.gene)
    assert gc_matrix_raw.data.shape[1] == len(gc_matrix_raw.cell)

    return gc_matrix_raw

def upload_dataset_atac_format(
    dataset, fdata, frow, fcol, database=DATABASE_ATAC):
    """
    Example:
        dataset = 'CEMBA_1B_180118'
        fdir = '/cndd/Public_Datasets/CEMBA/snATACSeq/Datasets/1B/CEMBA180118_1B_rep1/counts/genebody/'
        fdata = fdir + 'CEMBA180118_1B_rep1.genebody.count_matrix.npz'
        frow = fdir + 'CEMBA180118_1B_rep1.genebody.count_matrix.row.index' # gene
        fcol = fdir + 'CEMBA180118_1B_rep1.genebody.count_matrix.col.index' # col
        database = 'CEMBA_ATAC'

    """
    gc_matrix_raw = read_sparse_atac_matrix(fdata, frow, fcol, dataset)

    logging.info("Upload dataset: {}".format(dataset))
    upload_dataset_worker(dataset, gc_matrix_raw, database=DATABASE_ATAC)
    logging.info("Done Upload dataset: {}".format(dataset))

    return 

if __name__ == '__main__':
    log = create_logger()

    datasets_info = [
        # ('3C','CEMBA_3C_171206', 'CEMBA171206_3C_rep1'),
        # ('3C','CEMBA_3C_171207', 'CEMBA171207_3C_rep2'),
        # ('4B','CEMBA_4B_171212', 'CEMBA171212_4B_rep1'),
        # ('4B','CEMBA_4B_171213', 'CEMBA171213_4B_rep2'),
        # ('4B','CEMBA_4B_180104', 'CEMBA180104_4B_rep3'),
        # ('4D','CEMBA_4D_171214', 'CEMBA171214_4D_rep1'),
        # ('4D','CEMBA_4D_171219', 'CEMBA171219_4D_rep2'),
        # ('3F','CEMBA_3F_180105', 'CEMBA180105_3F_rep1'),
        # ('3F','CEMBA_3F_180109', 'CEMBA180109_3F_rep2'),
        ('4E','CEMBA_4E_180110', 'CEMBA180110_4E_rep1'),
        ('4E','CEMBA_4E_180111', 'CEMBA180111_4E_rep2'),
        # ('1B','CEMBA_1B_180118', 'CEMBA180118_1B_rep1'),
        # ('1B','CEMBA_1B_180119', 'CEMBA180119_1B_rep2'),
        # ('3A','CEMBA_3A_180129', 'CEMBA180129_3A_rep1'),
        # ('3A','CEMBA_3A_180130', 'CEMBA180130_3A_rep2'),
        # ('4A','CEMBA_4A_180205', 'CEMBA180205_4A_rep1'),
        # ('4A','CEMBA_4A_180206', 'CEMBA180206_4A_rep2'),
        # ('1C','CEMBA_1C_180208', 'CEMBA180208_1C_rep1'),
        # ('1C','CEMBA_1C_180212', 'CEMBA180212_1C_rep2'),
        # ('1B','CEMBA_1B_180213', 'CEMBA180213_1B_rep3'),
    ]

    for region, dataset, dataset_id in datasets_info:
        fdir = '/cndd/Public_Datasets/CEMBA/snATACSeq/Datasets/{}/{}/counts/genebody/'.format(region, dataset_id)
        fdata = fdir + '{}.genebody.count_matrix.npz'.format(dataset_id)
        frow = fdir + '{}.genebody.count_matrix.row.index'.format(dataset_id) # gene
        fcol = fdir + '{}.genebody.count_matrix.col.index'.format(dataset_id) # col

        upload_dataset_atac_format(dataset, fdata, frow, fcol, database=DATABASE_ATAC)












#!/usr/bin/env python3

import sqlalchemy as sa
from scipy import sparse
from collections import namedtuple

from __init__ import *
import snmcseq_utils
from snmcseq_utils import create_logger
from snmcseq_utils import cd 
from snmcseq_utils import sparse_logcpm 
from snmcseq_utils import sparse_logtpm
from CEMBA_update_mysql import insert_into
from CEMBA_update_mysql import connect_sql 
from CEMBA_update_mysql import gene_id_to_table_name 
from CEMBA_update_mysql import upload_to_datasets



# GC_matrix = namedtuple('GC_matrix', ['gene', 'cell', 'data'])
# def sparse_logcpm(gc_matrix):
#     """
#     """
#     lib_size_inv = sparse.diags(np.ravel(1.0/gc_matrix.data.sum(axis=0)))
#     logcpm = (gc_matrix.data).dot(lib_size_inv*1e6).tocoo()
#     logcpm.data = np.log10(logcpm.data + 1)

#     gc_logcpm = GC_matrix(
#         gc_matrix.gene, 
#         gc_matrix.cell, 
#         logcpm,
#     )
    
#     return gc_logcpm

# def sparse_logtpm(gc_matrix, gene_lengths):
#     """
#     gene_lengths: array like 
    
#     """
#     gene_lengths = np.array(gene_lengths)
#     gene_length_inv = sparse.diags(np.ravel(1.0/gene_lengths))
#     tmp = (gene_length_inv).dot(gc_matrix.data).tocoo()
#     lib_size_inv = sparse.diags(np.ravel(1.0/tmp.sum(axis=0)))
    
#     logtpm = tmp.dot(lib_size_inv*1e6).tocoo()
#     logtpm.data = np.log10(logtpm.data + 1)

#     gc_logtpm = GC_matrix(
#         gc_matrix.gene, 
#         gc_matrix.cell, 
#         logtpm,
#     )
    
#     return gc_logtpm

def upload_dataset_worker(dataset, gc_logtpm, database=DATABASE_RNA, 
    sex='N', brain_region='brain', target_region=None):
    """Given a gene-by-cell count matrix, calculate log10(CPM+1), upload it to CEMBA_RNA
    upload includes:
    - datasets table
    - cells table
    - gene_* tables
    """
    engine = connect_sql(database)

    # # logtpm
    # gc_logtpm = sparse_logtpm(gc_matrix_raw, gene_length)

    # upload_to_datasets
    upload_to_datasets(dataset, database=database, strict=False, 
        convention='CEMBA_RNA',
        sex=sex, 
        brain_region=brain_region, 
        target_region=target_region,
        )

    # upload to cells
    df_cells = pd.DataFrame()
    df_cells['cell_name'] = gc_logtpm.cell 
    df_cells['dataset'] = dataset 
    insert_into(engine, 'cells', df_cells, ignore=True, verbose='True')

    # upload to genes
    sql = """SELECT * FROM cells"""
    df_cells = pd.read_sql(sql, engine)
    metadata = sa.MetaData(engine)

    for i, gene in enumerate(gc_logtpm.gene):
        # if i > 10:
        #     break
        
        data_gene = gc_logtpm.data.getrow(i).tocoo()
        gene_table_name = gene_id_to_table_name(gene)

        if (i%100 == 0):
            logging.info("Progress on genes: {} {}".format(i+1, gene_table_name))
        
        df_gene = pd.DataFrame()
        df_gene['cell_name'] = [gc_logtpm.cell[col] for col in data_gene.col]
        df_gene['normalized_counts'] = data_gene.data
        
        df_gene = pd.merge(df_gene, df_cells, on='cell_name')
        df_gene = df_gene[['cell_id', 'normalized_counts']]

        # print(gene, df_gene)
        if not df_gene.empty:
            insert_into(engine, gene_table_name, df_gene, ignore=True, verbose=False)

    return

def upload_dataset_rna_format(dataset, filename, index_col, sep, gene_axis, 
    sex='N', brain_region='brain', target_region=None,
    cell_suffix=False,
    database=DATABASE_RNA, 
    path_genebody_annotation=PATH_GENEBODY_ANNOTATION):
    """Example:
        dataset = 'CEMBA_1B_180118'
        filename = ''
        sep = ',' or '\t'
        index_col = 'sample_id'
        gene_axis = 0 or 1
        database = 'CEMBA_ATAC'
        path

    """
    logging.info("To upload dataset: {}; sex: {}; brain_region: {}; target_region: {}; to {}\n"
            .format(dataset, sex, brain_region, target_region, database))

    # gene annotation
    df_genes = pd.read_table(path_genebody_annotation, index_col='gene_id')
    df_genes_v2 = df_genes.reset_index().groupby('gene_name').first()

    # input
    logging.info("Reading file from {}".format(filename))
    df = pd.read_csv(filename, sep=sep, index_col=index_col)
    logging.info("Done reading file")
    if gene_axis == 0:
        pass
    elif gene_axis == 1:
        df = df.T
    gene_names = df.index.values

    if cell_suffix:
        cells = ['{}_{}'.format(cell, dataset) for cell in df.columns.values]
    else:
        cells = df.columns.values

    # get non-duplicated gene list
    gene_ids = [snmcseq_utils.gene_name_to_id(name, df_genes_v2) for name in gene_names]
    good_idx = [i for i in range(len(gene_ids)) if gene_ids[i]]
    gene_names = np.array(gene_names)[good_idx]
    gene_ids = np.array(gene_ids)[good_idx]

    # truncate input
    df = df.loc[gene_names, :]
    df.index = gene_ids

    # structured 
    gc_counts = GC_matrix(
        gene=df.index, 
        cell=df.columns,
        data=sparse.coo_matrix(df.values),
        )
    df = None

    # tpm
    gene_lengths = (df_genes['end'] - df_genes['start']).loc[gc_counts.gene].values
    gc_logtpm = sparse_logtpm(gc_counts, gene_lengths)

    # upload  
    logging.info("Upload dataset: {}".format(dataset))
    upload_dataset_worker(dataset, gc_logtpm, database=DATABASE_RNA, 
                    sex=sex, brain_region=brain_region, target_region=target_region)
    logging.info("Done Upload dataset: {}".format(dataset))

    return 

if __name__ == '__main__':
    log = create_logger()

    datasets_info = [
        {'dataset': 'Zeng_SMARTer_nuclei_MOp', 
        'filename': '/cndd/Public_Datasets/CEMBA/BICCN_minibrain_data/Zeng/SMARTer_nuclei_MOp/exon_counts.csv.gz', 
        'index_col': 'sample_id', 
        'sep': ',', 
        'gene_axis': 1, 
        'sex': 'N', 
        'brain_region': 'MOp', 
        'target_region': None,
        'cell_suffix': False,
        },
    ]

    for dataset_info in datasets_info:
        upload_dataset_rna_format(
                dataset_info['dataset'], 
                dataset_info['filename'], 
                dataset_info['index_col'], 
                dataset_info['sep'], 
                gene_axis=dataset_info['gene_axis'], 
                sex=dataset_info['sex'], 
                brain_region=dataset_info['brain_region'], 
                target_region=dataset_info['target_region'],
                cell_suffix=dataset_info['cell_suffix'],
                database=DATABASE_RNA, 
                path_genebody_annotation=PATH_GENEBODY_ANNOTATION,
                )









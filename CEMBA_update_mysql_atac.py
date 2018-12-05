#!/usr/bin/env python3
"""Future Fangming: code to upload smoothed counts is un-tested.
"""

from __init__ import *
import sqlalchemy as sa
from scipy import sparse
from collections import namedtuple

import snmcseq_utils
import mnni_utils

from snmcseq_utils import create_logger
from snmcseq_utils import cd 
from snmcseq_utils import sparse_logcpm 
from snmcseq_utils import sparse_logtpm
from CEMBA_update_mysql import insert_into
from CEMBA_update_mysql import insert_into_worker
from CEMBA_update_mysql import connect_sql 
from CEMBA_update_mysql import gene_id_to_table_name 
from CEMBA_update_mysql import upload_to_datasets


def upload_atac_worker(dataset, dataset_cemba, database=DATABASE_ATAC, 
                    ens_id=None, 
                    ensemble_table=None, 
                    gc_matrix=None,
                    table_cells=True, table_datasets=True, 
                    table_ensembles=True, table_ens=True, 
                    table_genes=True,
                    ):
    """
    """
    # assert ensemble_table['cell_name'].tolist() == gc_matrix.cell.tolist()

    engine = connect_sql(database)

    # upload_to_datasets
    if table_datasets:
        upload_to_datasets(dataset, database=database, strict=False)

    # upload to cells
    if table_cells:
        df_cells = pd.DataFrame()
        df_cells['cell_name'] = ensemble_table['cell_name'] 
        df_cells['dataset'] = dataset 
        # if a cell barcode doesn't exist, it will be uploaded
        insert_into(engine, 'cells', df_cells, ignore=True, verbose=True) 

    # get cell_ids
    sql = """SELECT cell_id, cell_name FROM cells WHERE dataset = '{}'""".format(dataset)
    df_cells = pd.read_sql(sql, engine)
    # print(df_cells.head())

    # update to ensembles
    if table_ensembles:
        table_name = 'ensembles'
        ens_name = dataset
        ens_datasets = dataset
        snmc_ens_id = ens_id
        dict_list = [{'ensemble_id': ens_id, 
                      'ensemble_name': ens_name, 
                      'public_access': False, 
                      'datasets': ens_datasets, 
                      'snmc_ensemble_id': snmc_ens_id, 
                     }]
        # print(dict_list)
        insert_into_worker(engine, table_name, dict_list, ignore=False)

    # create and upload to Ens table
    if table_ens:
        engine.execute("DROP TABLE IF EXISTS Ens{}".format(ens_id))
        engine.execute("CREATE TABLE Ens{} LIKE Ens".format(ens_id))

        table_name = 'Ens{}'.format(ens_id)
        ensemble_table = pd.merge(ensemble_table, df_cells, on='cell_name').drop('cell_name', axis=1)
        # print(ensemble_table.head())
        insert_into(engine, table_name, ensemble_table, ignore=False, verbose=True)

    # upload to gene tables
    if table_genes:
        for i, gene in enumerate(gc_matrix.gene):
            gene_table_name = gene_id_to_table_name(gene)

            data_gene = gc_matrix.data.getrow(i).tocoo() # no need to upload 0 values

            if (i%100 == 0):
                logging.info("Progress on genes: {} {}".format(i+1, gene_table_name))
            
            df_gene = pd.DataFrame()
            df_gene['cell_name'] = [gc_matrix.cell[col] for col in data_gene.col]
            df_gene['normalized_counts_renlab'] = data_gene.data

            df_gene = pd.merge(df_gene, df_cells, on='cell_name')
            df_gene = df_gene[['cell_id', 'normalized_counts_renlab']]

            # print(gene, df_gene.head())
            if not df_gene.empty:
                insert_into(engine, gene_table_name, df_gene, ignore=True, verbose=False)
            # break

    return 


def upload_dataset_worker(dataset, gc_normalized, gc_smoothed_normalized, database=DATABASE_ATAC):
    """Given gene-by-cell count matrices (smoothed and not smoothed), upload it to CEMBA_ATAC
    """
    assert list(gc_normalized.gene) == list(gc_smoothed_normalized.gene)
    assert list(gc_normalized.cell) == list(gc_smoothed_normalized.cell)

    engine = connect_sql(database)

    # upload_to_datasets
    upload_to_datasets(dataset, database=database, strict=False)

    # upload to cells
    df_cells = pd.DataFrame()
    df_cells['cell_name'] = gc_normalized.cell 
    df_cells['dataset'] = dataset 
    insert_into(engine, 'cells', df_cells, ignore=True, verbose='True')

    # upload to genes
    sql = """SELECT * FROM cells"""
    df_cells = pd.read_sql(sql, engine)
    metadata = sa.MetaData(engine)

    for i, gene in enumerate(gc_normalized.gene):
        gene_table_name = gene_id_to_table_name(gene)

        data_gene = gc_normalized.data.getrow(i).tocoo()
        data_gene_smoothed = gc_smoothed_normalized.data.getrow(i).tocoo()

        if (i%100 == 0):
            logging.info("Progress on genes: {} {}".format(i+1, gene_table_name))
        
        df_gene = pd.DataFrame()
        df_gene['cell_name'] = [gc_normalized.cell[col] for col in data_gene.col]
        df_gene['normalized_counts'] = data_gene.data
        
        df_gene_smoothed = pd.DataFrame()
        df_gene_smoothed['cell_name'] = [gc_smoothed_normalized.cell[col] for col in data_gene_smoothed.col]
        df_gene_smoothed['smoothed_normalized_counts'] = data_gene_smoothed.data

        df_gene = pd.merge(df_gene, df_gene_smoothed, on='cell_name', how='outer')

        df_gene = pd.merge(df_gene, df_cells, on='cell_name')
        df_gene = df_gene[['cell_id', 'normalized_counts', 'smoothed_normalized_counts']]

        # print(gene, df_gene)
        if not df_gene.empty:
            insert_into(engine, gene_table_name, df_gene, ignore=True, verbose=False)
    return

def read_sparse_atac_matrix(fdata, frow, fcol, dataset='', add_dataset_as_suffix=True):
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

    if add_dataset_as_suffix:
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
    dataset, fdata, frow, fcol, fn_smooth_prefix=None, database=DATABASE_ATAC):
    """
    Example:
        dataset = 'CEMBA_1B_180118'
        fdir = '/cndd/Public_Datasets/CEMBA/snATACSeq/Datasets/1B/CEMBA180118_1B_rep1/counts/genebody/'
        fdata = fdir + 'CEMBA180118_1B_rep1.genebody.count_matrix.npz'
        frow = fdir + 'CEMBA180118_1B_rep1.genebody.count_matrix.row.index' # gene
        fcol = fdir + 'CEMBA180118_1B_rep1.genebody.count_matrix.col.index' # col
        database = 'CEMBA_ATAC'

    """
    # read in raw
    logging.info("Reading files ({})".format(fdata))
    gc_matrix_raw = read_sparse_atac_matrix(fdata, frow, fcol, dataset, add_dataset_as_suffix=True)

    # smooth and save
    logging.info("Smooth ATAC counts and save it to {}...".format(fn_prefix))
    df_genes = pd.read_table(PATH_GENEBODY_ANNOTATION, index_col='gene_id')
    gene_lengths = (df_genes['end'] - df_genes['start']).loc[gc_raw.gene]
    gc_logtpm = snmcseq_utils.sparse_logtpm(gc_raw, gene_lengths)

    p = 0.75 
    counts_atac_smoothed, M_aa = mnni_utils.smooth_in_modality(
                pd.DataFrame(gc_raw.data.todense(), index=gc_raw.gene, columns=gc_raw.cell),
                pd.DataFrame(gc_logtpm.data.todense(), index=gc_logtpm.gene, columns=gc_logtpm.cell),
                k=30, ka=5, npc=50, p=p)

    gc_counts_smoothed = GC_matrix(
                counts_atac_smoothed.index.values, 
                counts_atac_smoothed.columns.values, 
                sparse.coo_matrix(counts_atac_smoothed.values), 
                )

    if not fn_smooth_prefix:
        fn_prefix = '/cndd/Public_Datasets/CEMBA/snATACSeq/Datasets/{0}/{1}/counts_smoothed/genebody/{1}.genebody_smoothed'.format(region, dataset_id)
    save_gc_matrix(gc_counts_smoothed, fn_smooth_prefix)

    # normalize smoothed 
    gene_lengths = (df_genes['end'] - df_genes['start']).loc[gc_counts_smoothed.gene]
    gc_smoothed_logtpm = snmcseq_utils.sparse_logtpm(gc_counts_smoothed, gene_lengths)

    # upload
    logging.info("Upload dataset: {}".format(dataset))
    upload_dataset_worker(dataset, gc_logtpm, gc_smoothed_logtpm, database=DATABASE_ATAC)
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












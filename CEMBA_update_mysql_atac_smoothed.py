#!/usr/bin/env python3
"""Update/insert into smoothed_normalized_counts
"""

from __init__ import *
from scipy import sparse

import snmcseq_utils
import mnni_utils

import sqlalchemy as sa
from CEMBA_update_mysql import connect_sql
from CEMBA_update_mysql import gene_id_to_table_name 
from CEMBA_update_mysql import insert_into 
from CEMBA_update_mysql_atac import read_sparse_atac_matrix


def upload_smoothed_counts(dataset, gc_matrix, database=DATABASE_ATAC):
    """
    """
    engine = connect_sql(database)

    # upload to genes
    sql = """SELECT * FROM cells"""
    df_cells = pd.read_sql(sql, engine)
    metadata = sa.MetaData(engine)

    for i, gene in enumerate(gc_matrix.gene):
        
        data_gene = gc_matrix.data.getrow(i).tocoo()
        gene_table_name = gene_id_to_table_name(gene)

        if (i%1000 == 0):
            logging.info("Progress on genes: {} {}".format(i+1, gene_table_name))
        
        df_gene = pd.DataFrame()
        df_gene['cell_name'] = [gc_matrix.cell[col] for col in data_gene.col]
        df_gene['smoothed_normalized_counts'] = data_gene.data
        
        df_gene = pd.merge(df_gene, df_cells, on='cell_name')
        df_gene = df_gene[['cell_id', 'smoothed_normalized_counts']]
        
        sql = """SELECT {0}.cell_id FROM {0}
              JOIN cells ON cells.cell_id = {0}.cell_id
              WHERE cells.dataset = \'{1}\'""".format(gene_table_name, dataset)
        
        cell_ids_in_table = pd.read_sql(sql, engine)['cell_id'].values
        
        # if a cell is in the table, update the values; otherwise insert a row into the table 
        cond = df_gene['cell_id'].isin(cell_ids_in_table)
        
        # update table
        # df_gene[cond]
        gene_table = sa.Table(gene_table_name, metadata, autoload=True)
        for idx, row in df_gene[cond].iterrows():
            cell_id = int(row['cell_id'])
            smoothed_normalized_counts = float(row['smoothed_normalized_counts'])
            stmt = (gene_table.update()
                    .where(gene_table.c.cell_id==cell_id)
                    .values(smoothed_normalized_counts=smoothed_normalized_counts)
                   )
            stmt.execute()
        
        # insert into the table
        # df_gene[~cond]
        if not df_gene[~cond].empty:
            insert_into(engine, gene_table_name, df_gene[~cond], ignore=True, verbose=False)

    return


if __name__ == '__main__':

    log = snmcseq_utils.create_logger()
     
    # get atac file
    datasets_info = [
           #  ('3C','CEMBA_3C_171206', 'CEMBA171206_3C_rep1'),
            ('3C','CEMBA_3C_171207', 'CEMBA171207_3C_rep2'),
            ('4B','CEMBA_4B_171212', 'CEMBA171212_4B_rep1'),
            ('4B','CEMBA_4B_171213', 'CEMBA171213_4B_rep2'),
            ('4B','CEMBA_4B_180104', 'CEMBA180104_4B_rep3'),
            ('4D','CEMBA_4D_171214', 'CEMBA171214_4D_rep1'),
            ('4D','CEMBA_4D_171219', 'CEMBA171219_4D_rep2'),
            ('3F','CEMBA_3F_180105', 'CEMBA180105_3F_rep1'),
            ('3F','CEMBA_3F_180109', 'CEMBA180109_3F_rep2'),
            ('4E','CEMBA_4E_180110', 'CEMBA180110_4E_rep1'),
            ('4E','CEMBA_4E_180111', 'CEMBA180111_4E_rep2'),
            ('1B','CEMBA_1B_180118', 'CEMBA180118_1B_rep1'),
            ('1B','CEMBA_1B_180119', 'CEMBA180119_1B_rep2'),
            ('3A','CEMBA_3A_180129', 'CEMBA180129_3A_rep1'),
            ('3A','CEMBA_3A_180130', 'CEMBA180130_3A_rep2'),
            ('4A','CEMBA_4A_180205', 'CEMBA180205_4A_rep1'),
            ('4A','CEMBA_4A_180206', 'CEMBA180206_4A_rep2'),
            ('1C','CEMBA_1C_180208', 'CEMBA180208_1C_rep1'),
            ('1C','CEMBA_1C_180212', 'CEMBA180212_1C_rep2'),
            ('1B','CEMBA_1B_180213', 'CEMBA180213_1B_rep3'),
        ]

    for region, dataset, dataset_id in datasets_info:
        
        fp = '/cndd/Public_Datasets/CEMBA/snATACSeq/Datasets/{0}/{1}/counts_smoothed/genebody/{1}.genebody_smoothed'.format(region, dataset_id)
        fdata = fp + '.npz'
        frow = fp + '.row.index'
        fcol = fp + '.col.index'
        
        # try:
        logging.info("Reading files from dataset: {} \n({}*)".format(dataset, fp))
        gc_smoothed_counts = read_sparse_atac_matrix(fdata, frow, fcol, dataset, 
                                                     add_dataset_as_suffix=False)

        # logtpm normalize
        logging.info("logTPM normalize...")
        df_genes = pd.read_table(PATH_GENEBODY_ANNOTATION, index_col='gene_id')
        gene_lengths = (df_genes['end'] - df_genes['start']).loc[gc_smoothed_counts.gene]
        gc_smoothed_logtpm = snmcseq_utils.sparse_logtpm(gc_smoothed_counts, gene_lengths)
        
        logging.info("Uploading to mysql database: {}".format(DATABASE_ATAC))
        upload_smoothed_counts(dataset, gc_smoothed_logtpm, database=DATABASE_ATAC)
        # except:
        # logging.info("Problems with dataset {}, skipped!".format(dataset))
        



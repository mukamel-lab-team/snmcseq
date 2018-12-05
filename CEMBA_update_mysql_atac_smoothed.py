#!/usr/bin/env python3
"""Update/insert into smoothed_normalized_counts
"""

from __init__ import *
from scipy import sparse
import subprocess as sp 

import snmcseq_utils
import mnni_utils

import sqlalchemy as sa
from sqlalchemy import text 
from CEMBA_update_mysql import connect_sql
from CEMBA_update_mysql import gene_id_to_table_name 
from CEMBA_update_mysql import insert_into 
from CEMBA_update_mysql_atac import read_sparse_atac_matrix


def upload_smoothed_counts(dataset, gc_matrix, feature_type='smoothed_normalized_counts', database=DATABASE_ATAC):
    """
    """
    engine = connect_sql(database)

    # upload to genes
    sql = """SELECT * FROM cells WHERE dataset = '{}'""".format(dataset)
    df_cells = pd.read_sql(sql, engine)
    metadata = sa.MetaData(engine)
    csr_data = gc_matrix.data.tocsr() # important for slicing by rows later 

    for i, gene in enumerate(gc_matrix.gene):
        
        ti = time.time()
        # prepare the table to be uploaded
        data_gene = csr_data.getrow(i).tocoo()
        gene_table_name = gene_id_to_table_name(gene)

        if (i%1000 == 0):
            logging.info("Progress on genes: {} {}".format(i+1, gene_table_name))
        
        df_gene = pd.DataFrame()
        df_gene['cell_name'] = [gc_matrix.cell[col] for col in data_gene.col]
        df_gene[feature_type] = data_gene.data
        
        df_gene = pd.merge(df_gene, df_cells, on='cell_name')
        df_gene = df_gene[['cell_id', feature_type]]
        

        if not df_gene.empty: # otherwise just skip; the check is necessary to prevent error
            tmp_table = 'tmp_gene'
            full_table = gene_table_name

            # create a tmp_table
            sqls = [
                text('DROP TABLE IF EXISTS {}'.format(tmp_table)),
                text('CREATE TABLE {} LIKE templates.gene_table_atac_smoothed'.format(tmp_table)),
                text('ALTER TABLE {} CHANGE COLUMN `smoothed_normalized_counts` `{}` FLOAT'.format(tmp_table, feature_type)),
                ]
                # ]
            for sql in sqls:
                engine.execute(sql)

            # upload df_gene to the tmp_table
            insert_into(engine, tmp_table, df_gene, ignore=False, verbose=False)

            # insert tmp_table into the full gene table ON DUPLICATED KEY UPDATE
            sqls = [
                text('INSERT INTO {0} (cell_id, {2}) '
                     'SELECT cell_id, {2} from {1} '
                     'ON DUPLICATE KEY UPDATE {2} = VALUES({2})'
                     .format(full_table, tmp_table, feature_type)
                    ),
                text('DROP TABLE {}'.format(tmp_table)), 
                ] 
            for sql in sqls:
                engine.execute(sql)

            # done
    return

# deprecated
# def upload_smoothed_counts_old(dataset, gc_matrix, database=DATABASE_ATAC):
#     """
#     """
#     engine = connect_sql(database)

#     # upload to genes
#     sql = """SELECT * FROM cells WHERE dataset = '{}'""".format(dataset)
#     df_cells = pd.read_sql(sql, engine)
#     metadata = sa.MetaData(engine)
#     csr_data = gc_matrix.data.tocsr() # improve the speed later

#     for i, gene in enumerate(gc_matrix.gene):
        
#         ti = time.time()
#         # prepare the table to be uploaded
#         print('---'*5)
#         print(time.time()-ti)
#         data_gene = csr_data.getrow(i).tocoo()
#         print(time.time()-ti)
#         gene_table_name = gene_id_to_table_name(gene)
#         print(time.time()-ti)

#         if (i%1000 == 0):
#             logging.info("Progress on genes: {} {}".format(i+1, gene_table_name))
        
#         df_gene = pd.DataFrame()
#         df_gene['cell_name'] = [gc_matrix.cell[col] for col in data_gene.col]
#         df_gene['smoothed_normalized_counts'] = data_gene.data
        
#         df_gene = pd.merge(df_gene, df_cells, on='cell_name')
#         df_gene = df_gene[['cell_id', 'smoothed_normalized_counts']]
        
#         tmp_table = 'tmp_gene'
#         tmp_file = '/scratch/CEMBA_update_mysql_atac_smoothed_tmp_file.tsv'
#         full_table = gene_table_name

#         print(time.time()-ti)

#         # save the file locally
#         df_gene.to_csv(tmp_file, sep='\t', na_rep='NA', header=False, index=False)
#         print(time.time()-ti)

#         # create a tmp_table
#         sqls = [
#             text('DROP TABLE IF EXISTS {}'.format(tmp_table)),
#             text('CREATE TABLE {} LIKE templates.gene_table_atac_smoothed'.format(tmp_table)),
#             # text("LOAD DATA LOCAL INFILE '{}' INTO TABLE {}".format(tmp_file, tmp_table)), # this is not supported by pymysql
#             ]
#         for sql in sqls:
#             engine.execute(sql)

#         print(time.time()-ti)

#         # upload file to the tmp_table (not supported by pymysql; use bash as a work-around)
#         sql = "LOAD DATA LOCAL INFILE '{}' INTO TABLE {}".format(tmp_file, tmp_table)
#         cmd = "mysql -h {} -u {} -p{} {} -e \"{}\"".format(HOST, USER, PWD, DATABASE_ATAC, sql)
#         res = sp.run(cmd, shell=True)
#         if res.returncode != 0:
#           raise ValueError("Bash command returned with error!")

#         print(time.time()-ti)
#         # insert tmp_table into the full gene table ON DUPLICATED KEY UPDATE
#         sqls = [
#             text('INSERT INTO {} (cell_id, smoothed_normalized_counts) '
#                  'SELECT cell_id, smoothed_normalized_counts from {} '
#                  'ON DUPLICATE KEY UPDATE smoothed_normalized_counts = VALUES(smoothed_normalized_counts)'
#                  .format(full_table, tmp_table)
#                 ),
#             # text('DROP TABLE {}'.format(tmp_table))
#             ] 
#         for sql in sqls:
#             engine.execute(sql)
#         print(time.time()-ti)

#         # done
#     return


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
        
        try:
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

        except:
            logging.info("Problems with dataset {}, skipped!".format(dataset))
        



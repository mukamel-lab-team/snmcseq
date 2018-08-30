#!/usr/bin/env python3

import sqlalchemy as sa
from sqlalchemy import update
from sqlalchemy.dialects.mysql import INTEGER
from sqlalchemy.dialects.mysql import SMALLINT
from collections import OrderedDict
from itertools import product 


from __init__ import *
from snmcseq_utils import cd 
from snmcseq_utils import get_mouse_chromosomes
from snmcseq_utils import compute_global_mC 
from snmcseq_utils import create_logger
from CEMBA_run_tsne import run_tsne_CEMBA
from CEMBA_update_mysql import insert_into
from CEMBA_update_mysql import connect_sql 
from CEMBA_update_mysql import gene_id_to_table_name 

def gene_id_to_name(gene_id, df_genes):
    """df_genes
    """
    return df_genes.loc[gene_id, 'gene_name']

def gene_name_to_id(gene_name, df_genes):
    """df_genes
    """
    return df_genes[df_genes['gene_name'] == gene_name].index.values[0]


if __name__ == '__main__':
    log = create_logger()
    database = 'CEMBA_snATAC'
    engine = connect_sql(database)
    dataset = 'CEMBA_3C_171207'

    sql = """SELECT * FROM genes"""
    df_genes = pd.read_sql(sql, engine, index_col='gene_id')
    df_genes.head()

    with cd('/cndd/Public_Datasets/CEMBA/snATACSeq/MiniBrain_fromRongxin'):
        df_clst = pd.read_table('XW46_cluster.txt', header=None, names=['sample', 'cluster_ID'], index_col='sample')
        df_tsne = pd.read_table('XW46_tsne.txt', header=None, names=['sample', 'tsne_x', 'tsne_y'], index_col='sample')
        df = pd.read_table('XW46_gene.txt.gz') 

    df = df.T

    df_res = pd.merge(df_genes[['gene_name']], df, left_on='gene_name', right_index=True)
    df_res = df_res.drop('gene_name', axis=1)

    sql = """SELECT * FROM cells"""
    df_cells = pd.read_sql(sql, engine)

    for i, (gene_id, row) in enumerate(df_res.iterrows()):
        gene_table_name = gene_id_to_table_name(gene_id)
        if (i%1000==0):
            logging.info("Progress on genes: {} {}".format(i+1, gene_table_name))
        df_gene = row.to_frame('normalized_counts')
        df_gene.index = dataset + '_' + df_gene.index.values
        df_gene = pd.merge(df_gene, df_cells, left_index=True, right_on='cell_name', how='left')
        df_gene = df_gene[['cell_id', 'normalized_counts']]

        insert_into(engine, gene_table_name, df_gene, ignore=True, verbose=False)

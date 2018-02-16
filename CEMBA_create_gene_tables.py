#!/usr/bin/env python3

import pandas as pd 
import numpy as np
import sqlalchemy as sa
from sqlalchemy.dialects.mysql import INTEGER
from sqlalchemy.dialects.mysql import SMALLINT
import time
import os
import logging
import glob


from __init__ import *
from snmcseq_utils import compute_global_mC
from snmcseq_utils import create_logger

def gene_id_to_table_name(gene_id):
    """
    """
    table_name = 'gene_' + gene_id.replace('.', '_')
    return table_name


CELLS_TABLE_COLS = ['cell_id', 
                     'cell_name', 
                     'dataset', 
                     'cell_type',
                     'global_mCH', 
                     'global_mCG',
                     'global_mCA',
                     'global_mCCC', 
                     'estimated_mCH', 
                     'estimated_mCG',
                     'percent_genome_covered', 
                     'total_reads',
                     'mapped_reads', 
                     'mapping_rate', 
                     'nonclonal_reads', 
                     'percent_nonclonal_rate',
                     'filtered_reads',
                     'filtered_rate',
                     'lambda_mC']

def define_cells_table(metadata):
    """
    """
    table = sa.Table('cells', metadata,
                sa.Column('cell_id', sa.Integer, primary_key=True, autoincrement=True), # auto increment and not null are implicitly defined
                sa.Column('cell_name', sa.String(255), nullable=False, unique=True), # unique but not primary key
                sa.Column('dataset', sa.String(40), nullable=False), # not in mapping summary
                sa.Column('cell_type', sa.String(20), nullable=True), # not in mapping summary
                     
                sa.Column('global_mCH', sa.Float, nullable=True),
                sa.Column('global_mCG', sa.Float, nullable=True),
                sa.Column('global_mCA', sa.Float, nullable=True), # not in mapping summary
                sa.Column('global_mCCC', sa.Float, nullable=True), 
                sa.Column('estimated_mCH', sa.Float, nullable=True),
                sa.Column('estimated_mCG', sa.Float, nullable=True),
                sa.Column('percent_genome_covered', sa.Float, nullable=True),
                     
                sa.Column('total_reads', sa.Integer, nullable=True),
                sa.Column('mapped_reads', sa.Integer, nullable=True),
                sa.Column('mapping_rate', sa.Float, nullable=True),
                sa.Column('nonclonal_reads', sa.Integer, nullable=True),
                sa.Column('percent_nonclonal_rate', sa.Float, nullable=True),
                sa.Column('filtered_reads', sa.Integer, nullable=True),
                sa.Column('filtered_rate', sa.Float, nullable=True),
                sa.Column('lambda_mC', sa.Float, nullable=True),
                )
    return table


def define_gene_table(metadata, gene_id):
    """
    """
    table_name = gene_id_to_table_name(gene_id)
    table = sa.Table(table_name, metadata,
                          # NOT auto increment not null primary key
                        sa.Column('cell_id', sa.Integer, sa.ForeignKey('cells.cell_id'), primary_key=True, autoincrement=False), 
                        sa.Column('mCH', INTEGER(unsigned=True), nullable=True),
                        sa.Column('CH', INTEGER(unsigned=True), nullable=True),
                        sa.Column('mCG', INTEGER(unsigned=True), nullable=True),
                        sa.Column('CG', INTEGER(unsigned=True), nullable=True),
                        sa.Column('mCA', INTEGER(unsigned=True), nullable=True),
                        sa.Column('CA', INTEGER(unsigned=True), nullable=True),
    )
    return table

def define_ens_table(metadata, ens):
    """
    """
    table = sa.Table(ens, metadata,
                          # NOT auto increment not null primary key
                        sa.Column('cell_id', sa.Integer, sa.ForeignKey('cells.cell_id'), primary_key=True, autoincrement=False), 
                        sa.Column('cluster', SMALLINT(unsigned=True), nullable=True),
                        sa.Column('annotation', sa.String(20), nullable=True),
                        sa.Column('tsne_x_mCH', sa.Float, nullable=False),
                        sa.Column('tsne_y_mCH', sa.Float, nullable=False),
                        sa.Column('tsne_x_mCG', sa.Float, nullable=True),
                        sa.Column('tsne_y_mCG', sa.Float, nullable=True),
                        sa.Column('tsne_x_mCA', sa.Float, nullable=True),
                        sa.Column('tsne_y_mCA', sa.Float, nullable=True),
    )
    return table

def define_enss_table(metadata):
    """
    """
    table = sa.Table('ensembles', metadata,
                          # NOT auto increment not null primary key
                        sa.Column('ensemble_id', sa.Integer, primary_key=True, autoincrement=False), 
                        sa.Column('ensemble_name', sa.String(255), nullable=False),
                        sa.Column('public_access', sa.Boolean, nullable=False, default=False),
    )
    return table

def define_genes_table(metadata):
    """
    """
    table = sa.Table('genes', metadata,
                          # NOT auto increment not null primary key
                        sa.Column('gene_id', sa.String(50), primary_key=True, autoincrement=False), 
                        sa.Column('gene_name', sa.String(100), nullable=False), 
                        sa.Column('chr', sa.String(5), nullable=True), 
                        sa.Column('start', sa.Integer, nullable=True), 
                        sa.Column('end', sa.Integer, nullable=True), 
                        sa.Column('strand', sa.CHAR(1), nullable=True), 
                        sa.Column('gene_type', sa.String(100), nullable=True), 
    )
    return table

def insert_into_worker(engine, table_name, dict_list, ignore=False):
    """
    """
    metadata = sa.MetaData(engine)
    table = sa.Table(table_name, metadata, autoload=True, autoload_with=engine)
    if ignore: 
        inserter = table.insert().prefix_with('IGNORE')
    else:
        inserter = table.insert()
        
    return engine.execute(inserter, dict_list)

def insert_into(engine, table_name, df_sql, ignore=False, verbose=True):
    """
    """
    # change null to none
    df_sql = df_sql.where(pd.notnull(df_sql), None)
    if verbose:
        print(df_sql.head())
    # dict_list
    dict_list = list(df_sql.T.to_dict().values())
    
    return insert_into_worker(engine, table_name, dict_list, ignore=ignore)
    

log = create_logger()
# An engine connects to a mysql database
engine = sa.create_engine('mysql://f7xie:3405040212@localhost/CEMBA')


# create gene tables 
ti = time.time()

metadata = sa.MetaData(engine)
# load specfic tables needed to referenced as foreign key
sa.Table('cells', metadata, autoload=True)

# get all genes
sql = 'SELECT * FROM genes'
df_genes = pd.read_sql(sql, engine, index_col='gene_id')
# print(df_genes.shape)
# df_genes.head()
# df_genes.sort_index().head()

# create gene tables 
for i, gene_id in enumerate(df_genes.index):
    logging.info("Creating gene tables... ({}/{})".format(i+1, df_genes.index.shape[0]))
    gene_table = define_gene_table(metadata, gene_id)

metadata.create_all(engine)

tf = time.time()
logging.info("Done! Time spent: {} secends".format(tf - ti))
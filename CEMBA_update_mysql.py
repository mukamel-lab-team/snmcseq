#!/usr/bin/env python3

import pandas as pd 
import numpy as np
import sqlalchemy as sa
from sqlalchemy.dialects.mysql import INTEGER
import time
import os
import logging
import glob


from __init__ import *
from snmcseq_utils import create_logger



def gene_id_to_table_name(gene_id):
    """
    """
    table_name = 'gene_' + gene_id.replace('.', '_')
    return table_name

def connect_sql(database, user='f7xie', host='localhost', pwd='3405040212'):
	"""
	"""
	connect_string = 'mysql://{}:{}@{}/{}'.format(user, pwd, host, database)	
	return sa.create_engine(connect_string)


def define_cells_table(metadata):
    """
    """
    table = sa.Table('cells', metadata,
                sa.Column('cell_id', sa.Integer, primary_key=True, autoincrement=True), # auto increment and not null are implicitly defined
                sa.Column('cell_name', sa.String(255), nullable=False),
                sa.Column('dataset', sa.String(40), nullable=False),
                sa.Column('cell_type', sa.String(20), nullable=True),
                sa.Column('global_mCH', sa.Float, nullable=False),
                sa.Column('global_mCG', sa.Float, nullable=False),
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
                        sa.Column('cluster', sa.String(20), nullable=True),
                        sa.Column('annotation', sa.String(20), nullable=True),
                        sa.Column('tsne_x_mCH', sa.Float, nullable=False),
                        sa.Column('tsne_y_mCH', sa.Float, nullable=False),
                        sa.Column('tsne_x_mCG', sa.Float, nullable=True),
                        sa.Column('tsne_y_mCG', sa.Float, nullable=True),
                        sa.Column('tsne_x_mCA', sa.Float, nullable=True),
                        sa.Column('tsne_y_mCA', sa.Float, nullable=True),
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
 

def create_tables():
	"""
	Create cells table and gene tables

	This is not well written, dumped from jupyter
	"""
	ti = time.time()

	path_genebody = GENEBODY
	df_genes = pd.read_table(path_genebody, index_col='gene_id')

	metadata = sa.MetaData()
	# metadata.reflect(bind=engine)

	# define cells table
	cells_table = define_cells_table(metadata)

	# define gene tables 
	for i, gene_id in enumerate(df_genes.index):
	    logging.info("Creating gene tables... ({}/{})".format(i+1, df_genes.index.shape[0]))
	    gene_table = define_gene_table(metadata, gene_id)

	# define ensemble tables
	# ens = 'Ens1'
	# ens_table = define_ens_table(metadata, ens)

	metadata.create_all(engine)

	tf = time.time()
	print(tf - ti)

def metadata_to_cells_table(meta_path, dataset):
	"""
	"""
	df_meta = pd.read_table(meta_path)[['Sample', 'mCH/CH', 'mCG/CG']]
	df_cells = df_meta.copy()
	df_cells['cell_type'] = np.nan
	df_cells['dataset'] = dataset
	df_cells = df_cells[['Sample', 'dataset', 'cell_type', 'mCH/CH', 'mCG/CG']]
	df_cells.columns = [['cell_name', 'dataset', 'cell_type', 'global_mCH', 'global_mCG']]

	return df_cells


def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]

# def main():
# 	"""
# 	"""

if __name__ == '__main__':

	# main()	

	UPDATE_CELLS = False
	UPDATE_GENES = True

	log = create_logger()
	ti = time.time()

	database = 'CEMBA_test'
	dataset = 'CEMBA_3C_171206'	

	dataset_path = os.path.join(PATH_DATASETS, dataset)
	meta_path = os.path.join(dataset_path, 'mapping_summary_{}.tsv'.format(dataset)) 
	path_genebody = GENEBODY	
	genebody_paths = sorted(glob.glob(os.path.join(dataset_path, 'gene_level/genebody_*.tsv.bgz'.format(dataset))))	

	engine = connect_sql(database)

	if UPDATE_CELLS:
		# prepare updates for cells table
		df_cells = metadata_to_cells_table(meta_path, dataset)
		# update cells table
		insert_into(engine, 'cells', df_cells, ignore=True, verbose=False)


	if UPDATE_GENES:
		# update to gene tables
		logging.info("Begin updating gene tables...")

		# get current cells table from the database
		sql = """SELECT * FROM cells"""
		df_cells = pd.read_sql(sql, engine) 


		# getting genebody information (100 cells at a time) and upload to database
		n_max = 100
		for i, chunk_paths in enumerate(chunks(genebody_paths, n_max)):
			logging.info("Processing cells... ({}-{}/{})".format(i+1, i+n_max, len(genebody_paths)))

			dfs = []
			for genebody_path in chunk_paths:
				sample = os.path.basename(genebody_path)[len('genebody_'):-len('.tsv.bgz')]
				df = pd.read_table(genebody_path, index_col='gene_id', compression='gzip')
				df['sample'] = sample 
				dfs.append(df)

			dfs = pd.concat(dfs)

			for i, (gene_id, df_gene) in enumerate(dfs.groupby('gene_id')):
				if (i%1000==0):
					logging.info("Progress on genes: {}".format(i+1))
				df_gene = pd.merge(df_gene, df_cells, left_on='sample', right_on='cell_name', how='left')
				df_gene = df_gene[['cell_id', 'mCH', 'CH', 'mCG', 'CG', 'mCA', 'CA']]

				gene_table_name = gene_id_to_table_name(gene_id) 
				# updata df_gene to gene table
				insert_into(engine, gene_table_name, df_gene, ignore=True, verbose=False)


	tf = time.time()
	logging.info("Done updating dataset: {}".format(dataset))
	logging.info("Total time spent: {} sec".format(tf - ti))


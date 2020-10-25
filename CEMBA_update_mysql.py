#!/usr/bin/env python3

# import pandas as pd 
# import numpy as np
# import time
# import os
# import logging
# import glob

import argparse
import sqlalchemy as sa
from sqlalchemy.dialects.mysql import INTEGER
from sqlalchemy.dialects.mysql import SMALLINT
from collections import OrderedDict
from itertools import product

from __init__ import *
import snmcseq_utils
from snmcseq_utils import create_logger
from snmcseq_utils import get_chromosomes

from snmcseq_utils import compute_global_mC 
# import CEMBA_marker_genes_mysql

def gene_id_to_table_name(gene_id):
    """
    """
    table_name = 'gene_' + gene_id.replace('.', '_')
    return table_name

def connect_sql(database, user=USER, host=HOST, pwd=PWD):
    """
    """
    connect_string = 'mysql+pymysql://{}:{}@{}/{}'.format(user, pwd, host, database)    
    return sa.create_engine(connect_string)

def insert_into_worker(engine, table_name, dict_list, ignore=False):
    """
    """
    # convert float64 to float (update 09/25/2018 Fangming for pymysql)
    for dict_ in dict_list:
        for key in dict_.keys():
            if isinstance(dict_[key], np.float64):
                dict_[key] = float(dict_[key]) 
            elif isinstance(dict_[key], np.int64):
                dict_[key] = int(dict_[key]) 

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
        logging.info(df_sql.head())
    # dict_list
    dict_list = list(df_sql.T.to_dict().values())
    
    return insert_into_worker(engine, table_name, dict_list, ignore=ignore)
 

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]


def upload_to_datasets(dataset, database=DATABASE, strict=True, 
    convention=CONVENTION, 
    sex='N', brain_region='brain', target_region=None):
    """Add a entry (row) to datasets table
    """
    logging.info("Update dataset: {} to datasets table in {} database".format(dataset, database))
    import datetime

    datasets_table_content = []
    date = datetime.datetime.now().date()

    if not snmcseq_utils.isdataset(dataset):
        if strict:
            raise ValueError("May not be a valid dataset: {} (not found in snmCSeq datasets directory)".format(dataset))
        else:
            logging.warning("May not be a valid dataset: {} (not found in snmCSeq datasets directory)".format(dataset))
            
    if convention == 'CEMBA':
        if snmcseq_utils.isrs2(dataset): # rs2 dataset
            info = dataset.split('_')[2]
            injcode = info[0] 
            sex = info[1].upper()  
            slicecode = info[2:] 
            brain_region = snmcseq_utils.slicecode_to_region(slicecode)
            target_region = snmcseq_utils.injcode_to_region(injcode)
        else: # rs1 dataset
            sex = 'M'
            slicecode = dataset.split('_')[1] 
            brain_region = snmcseq_utils.slicecode_to_region(slicecode)
            target_region = None

    else:
        sex = sex
        brain_region = brain_region 
        target_region = target_region 

        
    datasets_table_content.append({
        'dataset': dataset,
        'date_online': date,
        'brain_region': brain_region,
        'target_region': target_region,
        'sex': sex, 
        })
    # insert
    engine = connect_sql(database)
    insert_into_worker(engine, 'datasets', datasets_table_content) 

    return

def upload_to_cells(dataset, database=DATABASE, 
                    update_mapping_summary=True):
    """
    """
    logging.info("Uploading cells from dataset {} to database {}".format(dataset, database))
    engine = connect_sql(database)

    dataset_path = os.path.join(PATH_DATASETS, dataset)
    meta_path = os.path.join(dataset_path, 'mapping_summary_{}.tsv'.format(dataset)) 

    df_meta = pd.read_table(meta_path)
    df_cells = df_meta.copy()
    df_cells.columns = rename_ms_cols(df_meta.columns)

    # cols not in mapping summary
    df_cells['cell_type'] = np.nan
    df_cells['dataset'] = dataset
    df_res = compute_global_mC(dataset)
    df_cells['global_mCA'] = df_res.loc[df_cells.cell_name, 'global_mCA'].values
    # add by Fangming 2/20/2019
    if 'RS2' not in dataset and 'SCI' not in dataset:
        # fix by Fangming 9/12/2019
        cell_name_structure = df_cells['cell_name'].apply(lambda x: len(x.split('_'))).values
        if np.all(cell_name_structure >= 4): # all in the old
            try: # test if it's in column 12
                df_cells['NeuN'] = (df_cells['cell_name'].apply(lambda x: int(x.split('_')[-3][1:])!=12)   
                                                         .replace(False, '-')
                                                         .replace(True, '+')
                                    )
            except:
                raise ValueError("Cell names don't match with the convention (XXX_A12_XXX_XXX)")

        elif np.all(cell_name_structure == 2): # all in the new
            try: # test if it's in column 12
                df_cells['NeuN'] = (df_cells['cell_name'].apply(lambda x: int(x.split('_')[0].split('-')[-1][1:])!=12)   
                                                         .replace(False, '-')
                                                         .replace(True, '+')
                                    )
            except:
                raise ValueError("Cell names don't match with the convention (XXX_A12_XXX_XXX)")
    else:
        df_cells['NeuN'] = '+'

    df_cells = df_cells[CELLS_TABLE_COLS[1:]]

    # insert into
    insert_into(engine, 'cells', df_cells, ignore=False, verbose='True')
    logging.info("Upload complete.")

    # add global mCA to mapping summary 
    if update_mapping_summary:
        df_meta['mCA/CA'] = df_res.loc[df_meta.Sample, 'global_mCA'].values
        df_meta.to_csv(meta_path, sep='\t', na_rep='NA', header=True, index=False)
    return  

def upload_to_genes(dataset, database=DATABASE, species=SPECIES):
    """upload to gene tables
    """

    logging.info("Uploading to gene tables. Dataset: {}, Database: {}".format(dataset, database))
    ti = time.time()

    engine = connect_sql(database)
    # get current cells table from the database
    sql = """SELECT * FROM cells WHERE dataset='{}'""".format(dataset)
    df_cells = pd.read_sql(sql, engine) 


    chromosomes = tuple(['chr'+i for i in get_chromosomes(species)])
    sql = """SELECT * FROM genes WHERE chr IN {}""".format(chromosomes)
    df_genes = pd.read_sql(sql, engine, index_col='gene_id') 

    # info ordered by genes
    dict_genes = OrderedDict() 
    for gene_id in df_genes.index:
        dict_genes[gene_id] = []
        
    # getting genebody information (100 cells at a time) and organize it by genes 

    # read in table by cell (100 at a time)
    logging.info("Reading genebody info from cells... ({})".format(dataset))
    n_max = 100
    genebody_paths = sorted(glob.glob(os.path.join(PATH_DATASETS, dataset, 'gene_level', 'genebody_*.tsv.bgz')))
    if not genebody_paths: 
        raise ValueError("Unable to get genebody files")

    for i, chunk_paths in enumerate(chunks(genebody_paths, n_max)):
        logging.info("Progress: {}-{}/{}".format(i*n_max+1, (i+1)*n_max, len(genebody_paths)))
        dfs = []
        for genebody_path in chunk_paths:
            sample = os.path.basename(genebody_path)[len('genebody_'):-len('.tsv.bgz')]
            df = pd.read_table(genebody_path, index_col='gene_id', compression='gzip')
            df['sample'] = sample 
            dfs.append(df)

        # concat them
        dfs = pd.concat(dfs)
        # split by gene
        for i, (gene_id, df_gene) in enumerate(dfs.groupby('gene_id')):
            dict_genes[gene_id].append(df_gene)
            
    # concat by genes 
    logging.info("Concatenating info by genes...")
    dict_genes_temp = OrderedDict()
    for i, (gene_id, dfs_gene) in enumerate(dict_genes.items()):
        dict_genes_temp[gene_id] = pd.concat(dfs_gene)

    # upload 
    logging.info('Uploading to mySQL database: {}'.format(database))
    for i, (gene_id, df_gene) in enumerate(dict_genes_temp.items()):
        if (i%1000==0):
            logging.info("Progress on genes: {}".format(i+1))
        # merge with cells table to get cell_id
        df_gene = pd.merge(df_gene, df_cells, left_on='sample', right_on='cell_name', how='left')
        df_gene = df_gene[['cell_id', 'mCH', 'CH', 'mCG', 'CG', 'mCA', 'CA']]

        gene_table_name = gene_id_to_table_name(gene_id) 
        # updata df_gene to gene table
        insert_into(engine, gene_table_name, df_gene, ignore=True, verbose=False)

    tf = time.time()
    logging.info("Total time spent: {} sec".format(tf - ti))
    return


def upload_to_enss(ens, ens_name, ens_datasets, database=DATABASE):
    """
    """
    engine = connect_sql(database)
    table_name = 'ensembles'

    ens_id = int(ens[len('Ens'):])
    dict_list = [{'ensemble_id': ens_id, 'ensemble_name': ens_name, 'public_access': False, 'datasets': ','.join(ens_datasets)}]
    insert_into_worker(engine, table_name, dict_list, ignore=False)

    return  

def define_ens_table(metadata, ens):
    """
    """
    # prototypes
#   sa.Column('tsne_x_mCH', sa.Float, nullable=True),
#   sa.Column('tsne_y_mCH', sa.Float, nullable=True),
#   sa.Column('cluster', SMALLINT(unsigned=True), nullable=True),
#   sa.Column('annotation', sa.String(20), nullable=True),

    # tsnes
    contexts = CONTEXTS + ['m'.join(comb_contexts) for comb_contexts in COMBINED_CONTEXTS_LIST]
    perps = PERPLEXITIES
    tsne_types = ['m{}_ndim2_perp{}'.format(context, p) for (context, p) in product(contexts, perps)]
    
    # prepare columns
    args_tsne = ([sa.Column('tsne_x_{}'.format(tsne_type), sa.Float, nullable=True) 
            for tsne_type in tsne_types] 
            + [sa.Column('tsne_y_{}'.format(tsne_type), sa.Float, nullable=True) 
            for tsne_type in tsne_types]
           )
    
    # clusters
    ks = K_NN
    cluster_types = ['m{}_lv_npc50_k{}'.format(context, k) for (context, k) in product(contexts, ks)]
    # prepare columns
    args_cluster = ([sa.Column('cluster_{}'.format(cluster_type), SMALLINT(unsigned=True), nullable=True) 
            for cluster_type in cluster_types] 
            + [sa.Column('annotation_{}'.format(cluster_type), sa.String(20), nullable=True) 
            for cluster_type in cluster_types]
           )
    
    # combine cols
    args = args_tsne + args_cluster
    # create the table
    table = sa.Table(ens, metadata,
                        # NOT auto increment not null primary key
                        sa.Column('cell_id', sa.Integer, sa.ForeignKey('cells.cell_id'), primary_key=True, autoincrement=False), 
                        *args
    )
    return table

def create_and_upload_to_ens(ens, database=DATABASE):
    """
    """
    logging.info("Creating {} table...".format(ens))
    engine = connect_sql(database)
    # create ens table
    metadata = sa.MetaData(engine)
    # load specfic tables needed to referenced as foreign key
    sa.Table('cells', metadata, autoload=True)
    # define ensemble table
    ens_table = define_ens_table(metadata, ens)
    metadata.create_all()
    logging.info("Done creating {} table!".format(ens))

    # upload to ens table
    logging.info("Uploading to {} table...".format(ens))
    ens_path = os.path.join(PATH_ENSEMBLES, ens)
    cluster_files = sorted(glob.glob(os.path.join(ens_path, 'cluster/cluster_*.tsv'))) 
    tsne_files = sorted(glob.glob(os.path.join(ens_path, 'tsne/tsne_*.tsv')))

    # get cluster and annotations
    dfs_cluster = []
    for cluster_file in cluster_files:
        # filename to cluster_type
        paras = os.path.basename(cluster_file).split('_')
        method, para1, para2, context = paras[1], paras[2], paras[3], paras[5]
        cluster_type = '{}_{}_{}_{}'.format(context, method, para1, para2) 
        # load data
        annot_file = cluster_file + '.annot'
        df_cluster = pd.read_table(cluster_file, index_col='sample')
        df_annot = pd.read_table(annot_file, index_col='cluster_ID')
        # merge 
        df_cluster = pd.merge(df_cluster, df_annot, left_on='cluster_ID', right_index=True)
        # organize format
        df_cluster.cluster_ID = [int(cluster_ID[len('cluster_'):])for cluster_ID in df_cluster.cluster_ID]
        df_cluster.columns = ['cluster_{}'.format(cluster_type), 
                             'annotation_{}'.format(cluster_type)]
        # append 
        dfs_cluster.append(df_cluster)
    # merge 
    df_clusters = pd.concat(dfs_cluster, axis=1)
        
    # get tsne
    dfs_tsne = []
    for tsne_file in tsne_files:
        # filename to cluster_type
        paras = os.path.basename(tsne_file).split('_')
        para1, para2, context = paras[1], paras[2], paras[5]
        tsne_type = '{}_{}_{}'.format(context, para1, para2) 
        # load data
        df_tsne = pd.read_table(tsne_file, index_col='sample')
        # organize format
        df_tsne.columns = [col+'_'+tsne_type for col in df_tsne.columns]
        # append
        dfs_tsne.append(df_tsne)
    df_tsnes = pd.concat(dfs_tsne, axis=1)
    df = pd.concat([df_tsnes, df_clusters], axis=1)


    # get cells
    sql = """SELECT * FROM cells"""
    df_cell = pd.read_sql(sql, engine)[['cell_id', 'cell_name']]

    # get cell_id
    df = pd.merge(df_cell, df, left_on='cell_name', right_index=True, how='right')
    df = df.drop('cell_name', axis=1)

    engine = connect_sql(database)
    table_name = ens 
    insert_into(engine, table_name, df, ignore=False, verbose=True)

    logging.info("Done uploading to {} table!".format(ens))

    return  

def upload_database_level(dataset, database=DATABASE):
    """
    """
    upload_to_datasets(dataset, database=database)
    upload_to_cells(dataset, database=database)
    upload_to_genes(dataset, database=database) 
    return

def upload_ensemble_level(ens, ens_name, ens_datasets, database=DATABASE):
    """
    """
    upload_to_enss(ens, ens_name, ens_datasets, database=database)
    create_and_upload_to_ens(ens, database=database)
    return



def create_parser():
    """
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-db", "--database", 
        required=True,
        help="database to upload to. Eg: CEMBA")
    parser.add_argument("-d", "--datasets", 
        nargs='*',
        help="List of dataset names. Eg: CEMBA_3C_171206")
    # parser.add_argument("-e", "--ensembles", 
    #     nargs='*',
    #     help="List of ensemble names. Eg: Ens1")

    return parser


if __name__ == '__main__':

    log = create_logger()
    ti = time.time()

    parser = create_parser()
    args = parser.parse_args()

    database = args.database 
    datasets = args.datasets
    # ensembles = args.ensembles

    if datasets:
        for dataset in datasets:
            upload_database_level(dataset, database=database)
    # if ensembles:
    #   for ens in ensembles:
    #       upload_ensemble_level(ens, ens_name, ens_datasets, database=database)

    tf = time.time()
    logging.info("Done updating dataset: {}".format(datasets))
    logging.info("Total time spent: {} sec".format(tf - ti))



#!/usr/bin/env python3
"""Upload analysis results from Rongxin Fang 12/4/2018
Upload his gene feature matrix
"""

from __init__ import *
import importlib
import pymysql

import snmcseq_utils
import CEMBA_update_mysql
import CEMBA_update_mysql_atac
import CEMBA_update_mysql_atac_smoothed
from scipy import sparse

logger = snmcseq_utils.create_logger()
logging.info('hello')

engine = CEMBA_update_mysql.connect_sql('CEMBA_ATAC')

# gene annotation
sql = 'SELECT * FROM genes'
df_genes = pd.read_sql(sql, engine).groupby('gene_name').first()
print(df_genes.head())

# Rongxin's analysis results
dirc = '/cndd/Public_Datasets/CEMBA/snATACSeq/fromRongxin_v3_20181130/analysis/*'
datasets_cemba = np.sort(np.unique([file.split('/')[-1].split('.')[0] for file in glob.glob(dirc)]))
datasets = np.array(['CEMBA_{}_{}'.format(dataset[-2:], dataset[5:5+6]) for dataset in datasets_cemba])
# datasets_cemba and datasets

# datasets available, ensembles available 
sql = "SELECT * FROM datasets"
datasets_available = pd.read_sql(sql, engine)['dataset'].values
sql = "SELECT * FROM ensembles JOIN datasets ON ensembles.datasets = datasets.dataset"
ensembles_available = pd.read_sql(sql, engine).set_index('dataset') #['dataset'].values

# datasets to upload
df_datasets = pd.DataFrame(datasets_cemba, index=datasets, columns=['dataset_cemba'])
df_datasets['exist'] = [1 if dataset in datasets_available else 0 for dataset in df_datasets.index]
df_datasets['ens_id'] = [ensembles_available.loc[dataset, 'ensemble_id'] 
                         if dataset in ensembles_available.index.values else -1 for dataset in df_datasets.index]
df_datasets = df_datasets[df_datasets['ens_id']!=-1]

# datasets
for dataset, dataset_cemba, *_, ens_id in df_datasets.itertuples():
    # upload one dataset
    print(dataset, dataset_cemba, ens_id)
    
    # prep info
    dirc = '/cndd/Public_Datasets/CEMBA/snATACSeq/fromRongxin_v3_20181130/analysis/'
    f1 = dirc + '{}.umap.txt'.format(dataset_cemba)
    df1 = pd.read_table(f1)
    f2 = dirc + '{}.cluster.txt'.format(dataset_cemba)
    df2 = pd.read_table(f2)
    df = pd.merge(df1, df2, left_index=True, right_index=True)
    df['cell_name'] = df['barcode'].apply(lambda x: '{}_{}'.format(x, dataset))

    ensemble_table = df.rename(columns={
              'umap1': 'tsne_x_ATAC', 
              'umap2': 'tsne_y_ATAC', 
              '0.5': 'cluster_ATAC', 
          })[['cell_name', 'tsne_x_ATAC', 'tsne_y_ATAC', 'cluster_ATAC']]
    
    f3 = dirc + '{}.gene.mat'.format(dataset_cemba)
    counts = pd.read_table(f3, index_col=0)
    counts = counts[np.intersect1d(counts.columns.values, df_genes.index.values)]
    counts.columns = [df_genes.loc[col, 'gene_id'] for col in counts.columns]
    counts.index = ['{}_{}'.format(cell, dataset) for cell in counts.index]

    counts = counts.loc[df['cell_name'].values]
    gc_matrix = GC_matrix(
        counts.columns.values, # gene 
        counts.index.values, # cell 
        sparse.coo_matrix(counts.values).T, # gene by cell
    )
    print('Begin upload')

    CEMBA_update_mysql_atac_smoothed.upload_smoothed_counts(
        dataset, gc_matrix, 
        feature_type='normalized_counts_renlab', 
        database=DATABASE_ATAC)

#     CEMBA_update_mysql_atac.upload_atac_worker(dataset, dataset_cemba, database=DATABASE_ATAC, 
#                     ens_id=ens_id, 
#                     ensemble_table=ensemble_table, 
#                     gc_matrix=None,
#                     table_datasets=False,
#                     table_cells=True,  
#                     table_ensembles=False, 
#                     table_ens=True, 
#                     table_genes=False,
#                     )

    # break

    
    



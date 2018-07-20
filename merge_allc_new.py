#!/usr/bin/env python3

from __init__ import *
import snmcseq_utils
import CEMBA_merge_allc
import CEMBA_update_mysql


log = snmcseq_utils.create_logger()

f = '/cndd/kkolodzi/CEMBA/striatum/genebody/clustering/clustering_integrate_all_v8-5.tsv'
ens = 'Ens51'
threshold = 2 
nprocs = 4

# get cluster assignment
df_clst = pd.read_table(f, index_col=0) 
df_clst['cluster_ID'] = [int(clst[len('cluster_'):]) for clst in df_clst['cluster_ID']]
df_clst['count'] = 1

mcc_cells_all = []
output_files_all = []
# specify output_files
for clst, df_sub in df_clst.groupby('cluster_ID'):
    cells = df_sub.index
    mcc_cells = [cell for cell in cells if cell.endswith('_indexed')]
    if len(mcc_cells) >= threshold:
        mcc_cells_all.append(mcc_cells)
        output_files_all.append(os.path.join(PATH_ENSEMBLES, ens, 
                                             'allc_merged', 'allc_multimodal_v1_clst{}.tsv'.format(clst))) 

# get allc file location of each cell
sql = """SELECT cell_name, dataset FROM {}
    JOIN cells ON {}.cell_id = cells.cell_id
""".format(ens, ens)
engine = CEMBA_update_mysql.connect_sql('CEMBA') 

df_info = pd.read_sql(sql, engine, index_col='cell_name')
df_info['path'] = [os.path.join(PATH_DATASETS, dataset, 'allc', 'allc_{}.tsv.bgz'.format(cell)) 
                   for cell, dataset in zip(df_info.index, df_info['dataset'])]

# get allc_files 
allc_files_all = []
for mcc_cells in mcc_cells_all:
    allc_files_all.append([df_info.loc[cell, 'path'] for cell in mcc_cells])

# do merging
CEMBA_merge_allc.merge_allc_parallel(allc_files_all, output_files_all, nprocs=nprocs)

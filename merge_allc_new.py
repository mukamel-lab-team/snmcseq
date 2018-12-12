#!/usr/bin/env python3

from __init__ import *
import snmcseq_utils
import CEMBA_merge_allc
import CEMBA_update_mysql


log = snmcseq_utils.create_logger()

f = '/cndd/fangming/human_mouse/results_keep/cluster_v2-181120.tsv'
f_meta = '/cndd/fangming/human_mouse/human_metadata_3human.tsv'
ens = 'Ens5'
nprocs = 4
threshold = 2 # number of cells included 

# get cluster assignment
df_clst = pd.read_table(f, index_col=0) 
df_clst['cluster_ID'] = [int(clst[len('cluster_'):]) for clst in df_clst['cluster_ID']]
df_clst['count'] = 1
# print(df_clst.shape)

# get cell info
df_meta = pd.read_table(f_meta, index_col=0) # get cells in the modality
# print(df_meta.shape)
df_clst =  pd.merge(df_clst, df_meta, left_index=True, right_index=True)
# print(df_clst.shape)

mcc_cells_all = []
output_files_all = []
# specify output_files
for clst, df_sub in df_clst.groupby('cluster_ID'):
    cells = df_sub.index
    mcc_cells = cells.values.tolist()
    # mcc_cells = [cell for cell in cells if cell.endswith('_indexed')]
    if len(mcc_cells) >= threshold:
        mcc_cells_all.append(mcc_cells)
        output_files_all.append(os.path.join('/cndd/Public_Datasets/human_snmcseq/Ensembles', ens, 
                                             'allc_merged', 'allc_human_mouse_v2-181120_clst{}_human.tsv'.format(clst))) 
# summary table
output_summary = os.path.join('/cndd/Public_Datasets/human_snmcseq/Ensembles', ens,
                              'allc_merged', 'summary_allc_human_mouse_v2-181120_human.tsv')
df_clst.to_csv(output_summary, sep='\t', na_rep='NA', header=True, index=True)

# get allc file location of each cell
sql = """SELECT cell_name, dataset FROM {}
    JOIN cells ON {}.cell_id = cells.cell_id
""".format(ens, ens)
engine = CEMBA_update_mysql.connect_sql('human_snmcseq') 

df_info = pd.read_sql(sql, engine, index_col='cell_name')
# print(df_info.shape)
# print(df_meta.shape)
# print(df_clst.shape)
df_info['path'] = [os.path.join('/cndd/Public_Datasets/human_snmcseq/Datasets', dataset, 'allc', 'allc_{}.tsv.bgz'.format(cell)) 
                   for cell, dataset in zip(df_info.index.values, df_info['dataset'])]

# get allc_files 
allc_files_all = []
num = 0
for mcc_cells in mcc_cells_all:
    tmp = [df_info.loc[cell, 'path'] for cell in mcc_cells 
                                    if cell in df_info.index.values.tolist()
                                    ]
    allc_files_all.append(tmp)
    num += len(tmp)

print(num)
print(len(allc_files_all), len(output_files_all))
print(allc_files_all[0][:3], output_files_all[:3])

# do merging

CEMBA_merge_allc.merge_allc_parallel(allc_files_all, output_files_all, nprocs=nprocs)


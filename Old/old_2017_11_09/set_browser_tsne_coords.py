#!/usr/bin/env python3
import pandas as pd

# on banjo
# input_fname = '/srv/scmdb_py_newdata/data/human_combined/tsne_points_ordered.csv'
# input_cluster_fname = ('/cndd/fangming/side_project/snmcseq_dev/preprocessed/' 
# 					+ 'human_hv1_hv2_CH_100000_clusters.tsv')
# output_fname = '/srv/scmdb_py_newdata/data/human_combined/tsne_points_ordered.csv'
# output_fname = './tsne_points_ordered_2.csv'

# df = pd.read_table(input_fname, header=0)
# df_cluster = pd.read_table(input_cluster_fname, header=0)
# df_merge = pd.merge(df, df_cluster, left_on='samp', right_on='sample', how='left', sort=True)

# print(df.shape)
# print(df_cluster.shape)
# print(df_merge.shape)

# df_merge['cluster_final'] = df_merge['cluster_ID'].values
# df_merge['cluster_name'] = df_merge['cluster_ID'].values
# df_merge['cluster_label'] = df_merge['cluster_ID'].values
# df_merge['cluster_ordered'] = [int(item.strip('cluster_'))
# 						for item in df_merge['cluster_ID'].tolist()]
# df_merge['biosample'] = [int(item.split('_')[0].strip('hv')) 
# 						for item in df_merge['samp'].tolist()]

# df_merge = df_merge[['samp', 'tsne_x', 'tsne_y', 
# 				'cluster_final', 'cluster_name', 'cluster_label', 'cluster_ordered',
# 				'cluster_ortholog', 'biosample']]
# df_merge.sort_values('cluster_ordered', inplace=True)

# df_merge.to_csv(output_fname, sep='\t', header=True, index=False, na_rep='NA')

# print(df_merge.head())



input_fname = '/srv/scmdb_py_newdata/data/human_combined/tsne_points_ordered.csv'
input_cluster_fname = ('/cndd/fangming/side_project/snmcseq_dev/preprocessed/' 
					+ 'human_hv1_hv2_CH_100000_out_v1_25_rotated.tsv')
output_fname = '/srv/scmdb_py_newdata/data/human_combined/tsne_points_ordered_v2.csv'
# output_fname = './tsne_points_ordered_2.csv'
df = pd.read_table(input_fname, header=0)
columns = df.columns.tolist()
print(df.shape)
for element in ['tsne_x', 'tsne_y']:
	columns.remove(element)
df = df[columns]
print(df.shape)
df_cluster = pd.read_table(input_cluster_fname, header=0)
df_merge = pd.merge(df, df_cluster, left_on='samp', right_on='cells', how='left', sort=True)

print(df.head())
print(df_cluster.head())

# do someting
df_merge['tsne_x'] = df_merge['tsne_x'].values
df_merge['tsne_y'] = df_merge['tsne_y'].values


# output
df_merge = df_merge[['samp', 'tsne_x', 'tsne_y', 
				'cluster_final', 'cluster_name', 'cluster_label', 'cluster_ordered',
				'cluster_ortholog', 'biosample']]
df_merge.sort_values('cluster_ordered', inplace=True)

print(df_merge.head())
df_merge.to_csv(output_fname, sep='\t', header=True, index=False, na_rep='NA')



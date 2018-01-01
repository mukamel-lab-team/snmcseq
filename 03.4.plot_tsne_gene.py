#!/usr/bin/env python3
"""
"""
import pandas as pd
import os

from snmcseq_utils import get_id_from_name 
from snmcseq_utils import plot_tsne_values 




def plot_tsne_gene_mcc(gene_name, 
					f_tsne, 
					context='CH', 
					normalized=False, 
					f_gene_id_to_names='/cndd/fangming/snmcseq_dev/data/references/gene_id_to_names.tsv',
					gene_mcc_dir=None,
					output=None,
					show=True,
					close=False):
	"""
	args: gene name, context, normalized/unnormalized -> a column in a gene_methylation file
			tsne coordinates
	"""
	if not gene_mcc_dir:
		if context == 'CH':
			gene_mcc_dir = '/cndd/fangming/snmcseq_dev/data/gene_level/genebody_mCH_human_combined_by_gene' 
		elif context == 'CG':
			gene_mcc_dir = '/cndd/fangming/snmcseq_dev/data/gene_level/genebody_mCG_human_combined_by_gene'

	# find gene mcc file
	gene_id = get_id_from_name(gene_name, f_gene_id_to_names=f_gene_id_to_names)
	f_gene_mcc = os.path.join(gene_mcc_dir, gene_id+'_m'+context+'.txt') 

	df_gene_mcc = pd.read_table(f_gene_mcc, header=None, index_col=1)

	df_gene_mcc.columns=['gene_id', 'mcc', 'nmcc']

	# get tsne coords
	df_tsne = pd.read_table(f_tsne, index_col='cells')

	# construct df_plot
	df_plot = pd.merge(df_tsne, df_gene_mcc, how='left', left_index=True, right_index=True) 
	if normalized:
		df_plot = df_plot[['tsne_x', 'tsne_y', 'nmcc']]
		title='Normalized genebody m%s: %s' % (context, gene_name)
		cbar_label='Normalized m%s' % context
		plot_tsne_values(df_plot, tx='tsne_x', ty='tsne_y', tc='nmcc', cbar_label=cbar_label,
				output=output, show=show, close=close, 
				t_xlim=None, t_ylim=None, title=title, figsize=(8,6))
	else:
		df_plot = df_plot[['tsne_x', 'tsne_y', 'mcc']]
		title='Genebody m%s: %s' % (context, gene_name)
		cbar_label='Normalized m%s' % context
		plot_tsne_values(df_plot, tx='tsne_x', ty='tsne_y', tc='mcc', cbar_label=cbar_label,
				output=output, show=show, close=close, 
				t_xlim=None, t_ylim=None, title=title, figsize=(8,6))





if __name__ == '__main__':


	# test example 
	# gene_name = 'GAD1'
	# f_tsne = './data/tsne/tsne_perp50_binc_mCHmCG_human_combined_100000_summary_normalized.tsv'
	# plot_tsne_gene_mcc(gene_name, 
	# 				f_tsne, 
	# 				context='CH', 
	# 				normalized=True, 
	# 				output='./test.pdf',
	# 				f_gene_id_to_names='/cndd/fangming/snmcseq_dev/data/references/gene_id_to_names.tsv',
	# 				gene_mcc_dir=None)


	gene_list = ['SATB2', 'TYRO3', 'ARPP21', 'SLC17A7', 'TBR1', 'CAMK2A', 'ITPKA', 
				'CUX1', 'CUX2', 'RORB', 'DEPTOR', 'VAT1L', 'SULF1', 'TLE4', 'FOXP2', 'GRIK3', 'BCL6', 
				'ERBB4','GAD1', 'SLC6A1', 
				'ADARB2', 'PROX1', 'SV2C', 
				'PVALB', 'SOX6', 'RELN', 'CACNA2D2', 'LHX6', 'GRIA1', 
				'MEF2C']

	f_tsne = './data/tsne/tsne_perp50_binc_mCHmCG_human_combined_100000_summary_normalized.tsv'
	output_dir='./results/tsne_marker_genes'
	context='CH'
	normalized=True
	for gene_name in gene_list:
		output = os.path.join(output_dir, 'tsne_marker_gene_%s_m%s_normalized.pdf' % (gene_name, context))
		plot_tsne_gene_mcc(gene_name, 
						f_tsne, 
						context=context, 
						normalized=normalized, 
						f_gene_id_to_names='/cndd/fangming/snmcseq_dev/data/references/gene_id_to_names.tsv',
						gene_mcc_dir=None,
						output=output,
						show=False,
						close=True)


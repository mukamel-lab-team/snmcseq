#!/usr/bin/env python3
import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np

BIN_SIZE = 100000
genome_path = '/cndd/fangming/iGenome/hg19/'
# genome_size_fname = '/cndd/fangming/iGenome/hg38/hg38.chrom.sizes'
output_path = '/cndd/fangming/side_project/snmcseq_dev/genome_stats'


def get_chrom_human(include_x=True, include_y=True):
	"""
	"""
	chromosomes = ['chr'+str(x) for x in range(1,23)]
	if include_x:
		chromosomes.append('chrX')
	if include_y:
		chromosomes.append('chrY')
	return chromosomes

def get_genome_size(genome_size_fname):
	"""
	"""
	srs_gsize = pd.read_table(genome_size_fname, header=None, index_col=0, squeeze=True) 
	return srs_gsize

def get_chrom_size_human(
	genome_size_fname='/cndd/fangming/iGenome/hg19/hg19.chrom.sizes'):	
	"""
	"""
	srs_gsize = get_genome_size(genome_size_fname=genome_size_fname)
	srs_gsize = srs_gsize.loc[get_chrom_human()]
	return srs_gsize




# main 

GEN_C_OCCURANCE = False 

# get total occurance of Cs for each bin
if GEN_C_OCCURANCE:
	# get human chroms
	chroms = get_chrom_human() 

	# get chrome sizes
	srs_gsize = get_chrom_size_human() 

	info = []
	for chrom in chroms:
		fname = os.path.join(genome_path, chrom+'.fa') 
		with open(fname, 'r') as seq_file:
			# remove the first line! >chrN
			seq_file.readline()
			# read in the rest of lines 
			seq = seq_file.read().replace('\n', '').upper()
			chrom_size = srs_gsize[chrom]
			# locally got seq length should be the same as chrom_size
			assert len(seq) == chrom_size # very important 

			bins = np.arange(0, chrom_size, BIN_SIZE)

			for bn in bins:
				if bn + BIN_SIZE > chrom_size: 
					break
				seq_bn = seq[bn:bn+BIN_SIZE]	
				info.append({'chrom': chrom,
							'bin': bn,
							'G_C': (seq_bn.count('C') + seq_bn.count('G')),
							})

	info = pd.DataFrame(info)
	info = info[['chrom', 'bin', 'G_C']]
	info.to_csv(os.path.join(genome_path, 'stats_binc.tsv'),
				sep='\t', na_rep='NA', header=True, index=False)
	print(info.head())


df_c = pd.read_table(os.path.join(genome_path, 'stats_binc.tsv'))

# get total number of covered Cs from allc file 

# ---

# get total number of covered Cs from 100kb bin matrix
input_coverage = '/cndd/fangming/side_project/snmcseq_dev/binc/human_combined_100kb_bins.tsv'
df_coverage = pd.read_table(input_coverage, header=0)

df_coverage['mean_hv1'] = df_coverage.filter(regex='^hv1_').mean(axis=1)
df_coverage['mean_hv2'] = df_coverage.filter(regex='^hv2_').mean(axis=1)

# # df_coverage['mean_tc'] = df_coverage.iloc[:, 2:].mean(axis=1)

# genome information
df_c = df_c[df_c.chrom != 'chrY']

df_plot = pd.merge(df_coverage, df_c, left_on=['chrom', 'bin'], right_on=['chrom', 'bin'])

df_plot['f_hv1'] = df_plot.mean_hv1/df_plot.G_C
df_plot['f_hv2'] = df_plot.mean_hv2/df_plot.G_C
df_plot['f_hv1'] = df_plot.f_hv1/df_plot.f_hv1.mean()
df_plot['f_hv2'] = df_plot.f_hv2/df_plot.f_hv2.mean()

print(df_plot.loc[((df_plot.f_hv1>2)|(df_plot.f_hv2>2)), 
	['chrom', 'bin', 'mean_hv1', 'mean_hv2']])


# not matched coordinates
# print(np.unique(df_coverage.chr.tolist()))
# print(np.unique(df_c.chrom.tolist()))
print(df_coverage.shape)
print(df_c.shape)
print(df_plot.shape)
# print(df_coverage.head())
# print(df_c.head())

# assert df_coverage.chr.tolist() == df_c.chrom.tolist()
# assert df_coverage.bin.tolist() == df_c.bin.tolist()
# a = df_coverage.mean_tc.values/df_c.C.values

# make a plot



## make plots
fig, axs = plt.subplots(2,1, figsize=(10,8))
ax = axs[0]
ax.plot(df_plot.f_hv1, label='human #1')
ax.plot(df_plot.f_hv2, label='human #2')
ax.set_title('Normalized coverage density of cytosine sites (whole genome)')
ax.set_ylabel('Normalized coverage density')
ax.legend()

ax = axs[1]
ax.plot(df_plot.f_hv1/df_plot.f_hv2, label='human #1/human #2')
# ax.set_title('Difference of normalized coverage density of cytosine sites (whole genome)')
ax.set_ylabel('Density ratio')
ax.legend()

fig.tight_layout()
output_fname = os.path.join(output_path, 'hs_totalc_coverage.pdf')
fig.savefig(output_fname)
print('Saved to %s' % output_fname)
plt.close(fig)
# plt.show()


for label, df_plot_sub in df_plot.groupby('chrom'): 
	fig, axs = plt.subplots(2, 1, figsize=(10,8))
	ax = axs[0]
	ax.plot(df_plot_sub.f_hv1, label='human #1')
	ax.plot(df_plot_sub.f_hv2, label='human #2')
	ax.legend()
	ax.set_title('Coverage density of cytosine sites: %s' %label)
	ax.set_ylabel('Normalized coverage density')

	ax = axs[1]
	ax.plot(df_plot_sub.f_hv1/df_plot_sub.f_hv2, label='human #1/human #2')
	# ax.set_title('Difference of normalized coverage density of cytosine sites: %s' %label)
	ax.set_ylabel('Density ratio')
	ax.legend()

	fig.tight_layout()
	output_fname = os.path.join(output_path, 'hs_totalc_coverage_%s.pdf' % label)
	fig.savefig(output_fname)
	print('Saved to %s' % output_fname)
	plt.close(fig)
	# plt.show()

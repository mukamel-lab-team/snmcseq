#!/usr/bin/env python3

from __init__ import *

f = '/cndd/Public_Datasets/CEMBA/snATACSeq/fromRongxin_v2_20180525/MOp_FPKM_gene_count_table.txt.gz'
df = pd.read_table(f)
print(df.shape)

dft = df.T
df = None
fo = '/cndd/Public_Datasets/CEMBA/snATACSeq/fromRongxin_v2_20180525/MOp_FPKM_gene_count_table_T.txt'
dft.to_csv(fo, sep='\t', na_rep='NA', header=False, index=True)
print('done')
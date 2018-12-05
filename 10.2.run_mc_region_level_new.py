#!/usr/bin/env python
"""
def run_mc_region_level(allc_files, output_files, 
	bed_file, 
	contexts=CONTEXTS,
	compress=True, 
	cap=2,
	nprocs=1):
"""

from __init__ import *
from CEMBA_run_mc_region_level import run_mc_region_level 

# allc_files = glob.glob('/cndd/Public_Datasets/CEMBA/snmCSeq/Ensembles/Ens100/allc_merged/allc_multimodal_v2_*.tsv.gz')[:2]
# bed_file = '/cndd/Public_Datasets/CEMBA/snmCSeq/Ensembles/Ens100/dmr/cgdmr_multimodal_v2_rms_results_collapsed.tsv.DMR_3dms.bed'
# contexts = ['CG'] # 'CH CG'

# clsts = [allc.split('/')[-1][len('allc_multimodal_v2_clst'):-len('.tsv.gz')] 
# 						for allc in allc_files]
# output_files = ['/cndd/Public_Datasets/CEMBA/snmCSeq/Ensembles/Ens100/dmr/cgcounts_cgdmr_multimodal_v2_clst{}.tsv'.format(clst) 
# 				for clst in clsts]

allc_files = glob.glob('/cndd/Public_Datasets/CEMBA/snmCSeq/Ensembles/Ens100/allc_merged/allc_multimodal_v2_*.tsv.gz')[:2]
bed_file = '/cndd/junhao/genomes/mm10/gencode.vM10.gene.annotation.bed'
contexts = ['CG'] # 'CH CG'
output_files = ['test_1.tsv', 'test_2.tsv']

cap = 0 # no counts cap
nprocs = 1
compress = True
bed_file_name_column = True 

print(allc_files, output_files, bed_file, contexts)

run_mc_region_level(allc_files, output_files, 
	bed_file, 
	bed_file_name_column=bed_file_name_column,
	contexts=contexts,
	compress=compress, 
	cap=cap,
	nprocs=nprocs)


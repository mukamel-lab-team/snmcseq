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
from natsort import natsorted

allc_files = natsorted(glob.glob(
	'/cndd/fangming/CEMBA/data/MOp_all/enhancers_may1/allc/allc_multimodal_v2_clst*.tsv.gz'
	))
bed_file = '/cndd/fangming/CEMBA/data/MOp_all/enhancers_may1/enhancers/allclusters_intersect_peaks_dmrs.bed.sort1'
contexts = ['CG'] # 'CH CG'
output_files = [('/cndd/fangming/CEMBA/data/MOp_all/enhancers_may1/enhancers_june19/counts/'
				+'mcg_'
				+allc.split('/')[-1][len('allc_'):-len('.tsv.gz')]+'.tsv') 
				for allc in allc_files]

cap = 0 # no counts cap
nprocs = 4
compress = True
bed_file_name_column = False

run_mc_region_level(allc_files, output_files, 
	bed_file, 
	bed_file_name_column=bed_file_name_column,
	contexts=contexts,
	compress=compress, 
	cap=cap,
	nprocs=nprocs)


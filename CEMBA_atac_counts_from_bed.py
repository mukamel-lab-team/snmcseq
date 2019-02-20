#!/usr/bin/env python

# Packages
import sys
sys.path.insert(0,'/cndd/widoyle/git_repos/CEMBA/mukamel_snATACSeq/')
import snatacseq_utils
import CEMBA_snATAC_preproc as preproc
import CEMBA_snATAC_qc as qcfunc
import CEMBA_snATAC_count as countfunc
import os
import glob
import subprocess
import shutil
import gzip

datasets_cemba = [
	# 'CEMBA171206_3C',
	'CEMBA171207_3C',
	'CEMBA171212_4B',
	'CEMBA171213_4B',
	]
	
for dataset_cemba in datasets_cemba:
	brain_slice = dataset_cemba.split('_')[1]
	readends = '/cndd/Public_Datasets/CEMBA/snATACSeq/Datasets/{1}/{0}/readends/{0}.readends.bed.gz'.format(dataset_cemba, brain_slice)
	barcodes = '/cndd/Public_Datasets/CEMBA/snATACSeq/Datasets/{1}/{0}/qc/{0}.qc_barcodes.txt'.format(dataset_cemba, brain_slice)
	bed = '/cndd/fangming/CEMBA/data/MOp_cicero/uniq_peaks_liftOver_nochr_mm10.bed'
	out = '/cndd/fangming/CEMBA/data/MOp_cicero/{}.cicero_peaks'.format(dataset_cemba)
	countfunc.counts_per_bed(readends,
	            barcodes,
	            bed,
	            out,
	            nproc=1,
	            combine=True)

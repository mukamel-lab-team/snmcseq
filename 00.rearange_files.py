#!/usr/bin/python3

# from __init__ import *
import subprocess as sp

PATH_ = '/cndd/Public_Datasets/CEMBA/snmCSeq/Datasets'

datasets = [
	'CEMBA_2A_180122',
	'CEMBA_2A_180123',
	'CEMBA_2C_180409',
	'CEMBA_2C_180410',
	'CEMBA_3B_180312',
	'CEMBA_3B_180501',
	'CEMBA_3D_180412',
	'CEMBA_3D_180416',
	'CEMBA_5B_180514',
	'CEMBA_5B_180529',
	'CEMBA_5D_180605',
	'CEMBA_5D_180612',
	'CEMBA_7B_180423',
	'CEMBA_7B_180424',
]

for dataset in datasets:
	a, b, c = dataset.split('_')
	dataset_cemba = 'CEMBA{}_{}'.format(c, b) 

	cmd = ("mkdir {2}/{0}/from_ecker_lab && " 
		"mv {2}/{3}/{1}/allc/*.gz {2}/{0}/from_ecker_lab"
		).format(dataset, dataset_cemba, PATH_, b)
	print(cmd)
	sp.run(cmd, shell=True)
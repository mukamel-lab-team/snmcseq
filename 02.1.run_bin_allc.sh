#!/bin/bash


input="/cndd/Public_Datasets/human_snmcseq/Datasets/MB_v1/allc \
		/cndd/Public_Datasets/human_snmcseq/Datasets/MB_EA/allc \
		/cndd/Public_Datasets/human_snmcseq/Datasets/MB_EB/allc  
		"

./CEMBA_run_bin_allc_files.py -f -i $input -n 8 


#!/bin/bash

# input="/cndd/Public_Datasets/CEMBA/snmCSeq/Datasets/CEMBA_4D_171214/allc \
#       /cndd/Public_Datasets/CEMBA/snmCSeq/Datasets/CEMBA_4D_171219/allc 
#      "

input="/cndd/Public_Datasets/human_snmcseq/Datasets/MB_v1/allc \
		/cndd/Public_Datasets/human_snmcseq/Datasets/MB_EA/allc \
		/cndd/Public_Datasets/human_snmcseq/Datasets/MB_EB/allc  
		"


./CEMBA_run_mc_gene_level.py -f -i $input -n 8 


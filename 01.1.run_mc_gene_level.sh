#!/bin/bash

input="/cndd/Public_Datasets/CEMBA/Datasets/CEMBA_3C_171206/allc \
	/cndd/Public_Datasets/CEMBA/Datasets/CEMBA_3C_171207/allc \
	/cndd/Public_Datasets/CEMBA/Datasets/CEMBA_4B_171212/allc \
	/cndd/Public_Datasets/CEMBA/Datasets/CEMBA_4B_171213/allc 
	"

./run_mc_gene_level.py -f -i $input -n 16

#!/bin/bash

samples="/cndd2/Public_Datasets/CEMBA_cndd2/snmCSeq/Datasets/to_process_gene_leftover.txt"

mapfile -t datasets < $samples
for data in ${datasets[@]}; do
	echo $data
	input="/cndd2/Public_Datasets/CEMBA_cndd2/snmCSeq/Datasets/$data/allc"
	echo $input
 	./CEMBA_run_mc_gene_level.py -f -chr -i $input -n 2 
done

#!/bin/bash

samples="/cndd2/Public_Datasets/CEMBA_cndd2/snmCSeq/Datasets/to_process_binc_leftover.txt"
mapfile -t datasets < $samples

for data in ${datasets[@]}; do
	echo $data
	input="/cndd2/Public_Datasets/CEMBA_cndd2/snmCSeq/Datasets/$data/allc"
	echo $input
	./CEMBA_run_bin_allc_files.py -f -i $input -n 4 
done


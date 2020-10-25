#!/bin/bash

database='CEMBA'

# space delimited string
samples="/cndd2/Public_Datasets/CEMBA_cndd2/snmCSeq/Datasets/to_process_upload2.txt"
mapfile -t datasets < $samples
datasets="${datasets[*]}"
echo $datasets

./CEMBA_update_mysql.py -db $database -d $datasets

#!/bin/bash


# input="/cndd/Public_Datasets/CEMBA/snmCSeq/Datasets/CEMBA_4C_180417/allc \
#        /cndd/Public_Datasets/CEMBA/snmCSeq/Datasets/CEMBA_4C_180419/allc
# "

input="/cndd/Public_Datasets/CEMBA/snmCSeq/Datasets/CEMBA_SCI_2017/allc \
"

./CEMBA_run_bin_allc_files.py -f -i $input -n 8 


#!/bin/bash

# input="/cndd/Public_Datasets/CEMBA/Datasets/CEMBA_3C_171206/allc \
# 	/cndd/Public_Datasets/CEMBA/Datasets/CEMBA_3C_171207/allc \
# 	/cndd/Public_Datasets/CEMBA/Datasets/CEMBA_4B_171212/allc \
# 	/cndd/Public_Datasets/CEMBA/Datasets/CEMBA_4B_171213/allc 
# 	"

# input="/cndd/Public_Datasets/CEMBA/snmCSeq/Datasets/CEMBA_4D_171214/allc \
#        /cndd/Public_Datasets/CEMBA/snmCSeq/Datasets/CEMBA_4D_171219/allc
# "

input="/cndd/Public_Datasets/CEMBA/snmCSeq/Datasets/CEMBA_3F_180109/allc"

./CEMBA_run_bin_allc_files.py -f -i $input -n 8 


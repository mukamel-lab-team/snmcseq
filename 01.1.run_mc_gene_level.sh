#!/bin/bash


input="/cndd/Public_Datasets/CEMBA/snmCSeq/Datasets/CEMBA_2A_180122/allc \
      /cndd/Public_Datasets/CEMBA/snmCSeq/Datasets/CEMBA_2A_180123/allc \
      /cndd/Public_Datasets/CEMBA/snmCSeq/Datasets/CEMBA_2C_180409/allc \
      /cndd/Public_Datasets/CEMBA/snmCSeq/Datasets/CEMBA_2C_180410/allc \
      /cndd/Public_Datasets/CEMBA/snmCSeq/Datasets/CEMBA_3B_180312/allc \
      /cndd/Public_Datasets/CEMBA/snmCSeq/Datasets/CEMBA_3B_180501/allc \
      /cndd/Public_Datasets/CEMBA/snmCSeq/Datasets/CEMBA_3D_180412/allc \
      /cndd/Public_Datasets/CEMBA/snmCSeq/Datasets/CEMBA_3D_180416/allc \
      /cndd/Public_Datasets/CEMBA/snmCSeq/Datasets/CEMBA_5B_180514/allc \
      /cndd/Public_Datasets/CEMBA/snmCSeq/Datasets/CEMBA_5B_180529/allc \
      /cndd/Public_Datasets/CEMBA/snmCSeq/Datasets/CEMBA_5D_180605/allc \
      /cndd/Public_Datasets/CEMBA/snmCSeq/Datasets/CEMBA_5D_180612/allc \
      /cndd/Public_Datasets/CEMBA/snmCSeq/Datasets/CEMBA_7B_180423/allc \
      /cndd/Public_Datasets/CEMBA/snmCSeq/Datasets/CEMBA_7B_180424/allc
     "

# input="/cndd/Public_Datasets/CEMBA/snmCSeq/Datasets/CEMBA_SCI_2017/allc \
# "

./CEMBA_run_mc_gene_level.py -f -i $input -n 4 


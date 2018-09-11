#!/bin/bash

# ens_id=5
# ens_name='CEMBA_3C'
# message="A combined ensemble with 2 datasets in CEMBA_3C"
# ens_sql="SELECT cell_name FROM cells WHERE dataset LIKE 'CEMBA_3C_%'"
# ./CEMBA_init_ensemble.py -ei $ens_id -en $ens_name -m $message \
# 						 --ensemble_sql $ens_sql

# ens_id=6
# ens_name='CEMBA_4B'
# message="A combined ensemble with 2 datasets in CEMBA_4B"
# # ens_datasets="CEMBA_4B_171212\
# # 			  CEMBA_4B_171213"
# ens_sql="SELECT cell_name FROM cells WHERE dataset LIKE 'CEMBA_4B_%'"
# # ens_cells=`cat ens_cells.txt`
# ./CEMBA_init_ensemble.py -ei $ens_id -en $ens_name -m $message \
# 						 --ensemble_sql $ens_sql

# ens_id=29
# ens_name='CEMBA_4D'
# message="A combined ensemble with 2 datasets in CEMBA 4D"
# ens_datasets="CEMBA_4D_171214\
# 			CEMBA_4D_171219"
			
# ./CEMBA_init_ensemble_v2.py -ei $ens_id -en $ens_name -m $message \
# 						 --ensemble_datasets $ens_datasets

ens_id=4
ens_name='HUMAN_MB_v1_EA_EB'
message="3 human snmcseq samples: MB_v1(published), MB_EA, MB_EB"
ens_datasets="MB_v1 \
MB_EA \
MB_EB"
			
./CEMBA_init_ensemble_v2.py -ei $ens_id -en $ens_name -m $message \
						 --ensemble_datasets $ens_datasets


# ens_id=31
# ens_name='CEMBA_RS2_MOp'
# message="This is an ensemble including all current RS2 datasets in MOp. (12 in total)"
# ens_datasets="CEMBA_RS2_Bm3C_rep1 \
# 	CEMBA_RS2_Bm3C_rep2 \
# 	CEMBA_RS2_Bm4B_rep1 \
# 	CEMBA_RS2_Bm4B_rep2 \
# 	CEMBA_RS2_Pf3C \
# 	CEMBA_RS2_Pf4B \
# 	CEMBA_RS2_Pm3C \
# 	CEMBA_RS2_Pm4B \
# 	CEMBA_RS2_Tf3C \
# 	CEMBA_RS2_Tf4B \
# 	CEMBA_RS2_Tm3C \
# 	CEMBA_RS2_Tm4B" 
# ./CEMBA_init_ensemble_v2.py -ei $ens_id -en $ens_name -m $message \
# 						 --ensemble_datasets $ens_datasets

# ens_id=32
# ens_name='CEMBA_MOp_with_RS2'
# message="This is an ensemble including all current datasets in MOp area, including RS2 datasets. (5+12=17 in total)"
# ens_datasets="CEMBA_3C_171206 \
# 	CEMBA_3C_171207 \
# 	CEMBA_4B_171212 \
# 	CEMBA_4B_171213 \
# 	CEMBA_4B_180104 \
# 	CEMBA_RS2_Bm3C_rep1 \
# 	CEMBA_RS2_Bm3C_rep2 \
# 	CEMBA_RS2_Bm4B_rep1 \
# 	CEMBA_RS2_Bm4B_rep2 \
# 	CEMBA_RS2_Pf3C \
# 	CEMBA_RS2_Pf4B \
# 	CEMBA_RS2_Pm3C \
# 	CEMBA_RS2_Pm4B \
# 	CEMBA_RS2_Tf3C \
# 	CEMBA_RS2_Tf4B \
# 	CEMBA_RS2_Tm3C \
# 	CEMBA_RS2_Tm4B" 
			
# ./CEMBA_init_ensemble_v2.py -ei $ens_id -en $ens_name -m $message \
# 						 --ensemble_datasets $ens_datasets

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


ens_id=117
ens_name='CEMBA_2A'
message="2A"
ens_datasets="CEMBA_2A_180122 \
CEMBA_2A_180123
"
			
./CEMBA_init_ensemble_v2.py -ei $ens_id -en $ens_name -m $message \
						 --ensemble_datasets $ens_datasets


ens_id=118
ens_name='CEMBA_2C'
message="2C"
ens_datasets="CEMBA_2C_180409 \
CEMBA_2C_180410
"
			
./CEMBA_init_ensemble_v2.py -ei $ens_id -en $ens_name -m $message \
						 --ensemble_datasets $ens_datasets


ens_id=119
ens_name='CEMBA_3B'
message="3B"
ens_datasets="CEMBA_3B_180312 \
CEMBA_3B_180501
"
			
./CEMBA_init_ensemble_v2.py -ei $ens_id -en $ens_name -m $message \
						 --ensemble_datasets $ens_datasets


ens_id=120
ens_name='CEMBA_3D'
message="3D"
ens_datasets="CEMBA_3D_180412 \
CEMBA_3D_180416
"
			
./CEMBA_init_ensemble_v2.py -ei $ens_id -en $ens_name -m $message \
						 --ensemble_datasets $ens_datasets

ens_id=121
ens_name='CEMBA_5B'
message="5B"
ens_datasets="CEMBA_5B_180514 \
CEMBA_5B_180529
"
			
./CEMBA_init_ensemble_v2.py -ei $ens_id -en $ens_name -m $message \
						 --ensemble_datasets $ens_datasets


ens_id=122
ens_name='CEMBA_5D'
message="5D"
ens_datasets="CEMBA_5D_180605 \
CEMBA_5D_180612
"
			
./CEMBA_init_ensemble_v2.py -ei $ens_id -en $ens_name -m $message \
						 --ensemble_datasets $ens_datasets


ens_id=123
ens_name='CEMBA_7B'
message="7B"
ens_datasets="CEMBA_7B_180423 \
CEMBA_7B_180424
"
			
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

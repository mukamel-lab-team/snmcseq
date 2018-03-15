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

# ens_id=7
# ens_name='CEMBA_MOp'
# message="A combined ensemble with 4 datasets in CEMBA 3C and 4B"
# ens_datasets="CEMBA_3C_171206\
# 			CEMBA_3C_171207\
# 			CEMBA_4B_171212\
# 			CEMBA_4B_171213"
# ./CEMBA_init_ensemble.py -ei $ens_id -en $ens_name -m $message \
# 						 --ensemble_datasets $ens_datasets

# ens_id=9
# ens_name='CEMBA_4B_v2'
# message="A combined ensemble with 3 datasets in CEMBA 4B"
# ens_datasets="CEMBA_4B_171212\
# 			CEMBA_4B_171213\
# 			CEMBA_4B_180104"

# ./CEMBA_init_ensemble_v2.py -ei $ens_id -en $ens_name -m $message \
# 						 --ensemble_datasets $ens_datasets


# ens_id=10
# ens_name='CEMBA_MOp_v2'
# message="A combined ensemble with 5 datasets in CEMBA 3C and 4B"
# ens_datasets="CEMBA_3C_171206\
# 			CEMBA_3C_171207\
# 			CEMBA_4B_171212\
# 			CEMBA_4B_171213\
# 			CEMBA_4B_180104"
			
# ./CEMBA_init_ensemble_v2.py -ei $ens_id -en $ens_name -m $message \
# 						 --ensemble_datasets $ens_datasets

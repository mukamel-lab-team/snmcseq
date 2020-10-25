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


# ens_id=154
# ens_name='CEMBA_11B'
# message="11B"
# ens_datasets="CEMBA_11B_190314 \
# CEMBA_11B_190325
# "
# 			
# ./CEMBA_init_ensemble_v2.py -ei $ens_id -en $ens_name -m $message \
#						 --ensemble_datasets $ens_datasets


declare -a ens_regions
ens_regions=(11F 12B 3F 4G 5A 5C 5E 5F 5G 5H 5J 6A 6B 6C 6D 8B 9B 9D 9H 9J)

ens_id=197 # last ens_id !!!
for ens_region in ${ens_regions[@]}
do
	ens_id=$((ens_id+1))
	ens_name="CEMBA_${ens_region}"
	message="A combined ensemble with 2 datasets in $ens_name"
	ens_sql="SELECT cell_name FROM cells WHERE dataset LIKE \"${ens_name}_%\""
	echo $ens_id $ens_name
	echo $message
	echo $ens_sql
	./CEMBA_init_ensemble_v2.py -ei $ens_id -en $ens_name -m $message \
				      --ensemble_sql $ens_sql
done



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

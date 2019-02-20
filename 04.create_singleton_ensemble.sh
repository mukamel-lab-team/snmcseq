#!/bin/bash

# dataset='MB_test'
# ens_id=1
# message="Singleton ensemble of $dataset dataset."
# ens_name=$dataset
# ./CEMBA_init_ensemble_v2.py --singleton -d $dataset -ei $ens_id -m $message -en $ens_name

# dataset='CEMBA_SCI_2017'
# ens_id=0
# message="Singleton ensemble of $dataset dataset."
# ens_name=$dataset
# ./CEMBA_init_ensemble_v2.py --singleton -d $dataset -ei $ens_id -m $message -en $ens_name

declare -A dataset_ids
dataset_ids['CEMBA_2A_180122']=77
dataset_ids['CEMBA_2A_180123']=78
dataset_ids['CEMBA_2C_180409']=112
dataset_ids['CEMBA_2C_180410']=113
dataset_ids['CEMBA_3B_180312']=83
dataset_ids['CEMBA_3B_180501']=114
dataset_ids['CEMBA_3D_180412']=84
dataset_ids['CEMBA_3D_180416']=85
dataset_ids['CEMBA_5B_180514']=96
dataset_ids['CEMBA_5B_180529']=97
dataset_ids['CEMBA_5D_180605']=102
dataset_ids['CEMBA_5D_180612']=103
dataset_ids['CEMBA_7B_180423']=115
dataset_ids['CEMBA_7B_180424']=116

for dataset in "${!dataset_ids[@]}"
do
	ens_id=${dataset_ids[$dataset]}
	# echo $dataset $ens_id
	message="Singleton ensemble of $dataset dataset."
	ens_name=$dataset
	./CEMBA_init_ensemble_v2.py --singleton -d $dataset -ei $ens_id -m $message -en $ens_name
done

# datasets="CEMBA_2A_180122 \
# 	CEMBA_2A_180123 \
# 	CEMBA_2C_180409 \
# 	CEMBA_2C_180410 \
# 	CEMBA_3B_180312 \
# 	CEMBA_3B_180501 \
# 	CEMBA_3D_180412 \
# 	CEMBA_3D_180416 \
# 	CEMBA_5B_180514 \
# 	CEMBA_5B_180529 \
# 	CEMBA_5D_180605 \
# 	CEMBA_5D_180612 \
# 	CEMBA_7B_180423 \
# 	CEMBA_7B_180424
# "

# ens_id=111 # LAST ensemble id!!! 
# for dataset in $datasets
# do
# 	ens_id=$((ens_id+1))
# 	# echo $dataset $ens_id
# 	message="Singleton ensemble of $dataset dataset."
# 	ens_name=$dataset
# 	./CEMBA_init_ensemble_v2.py --singleton -d $dataset -ei $ens_id -m $message -en $ens_name
# done

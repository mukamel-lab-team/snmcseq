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

# declare -A dataset_ids
# dataset_ids['CEMBA_2A_180122']=77
# dataset_ids['CEMBA_2A_180123']=78
# dataset_ids['CEMBA_2C_180409']=112
# dataset_ids['CEMBA_2C_180410']=113
# dataset_ids['CEMBA_3B_180312']=83
# dataset_ids['CEMBA_3B_180501']=114
# dataset_ids['CEMBA_3D_180412']=84
# dataset_ids['CEMBA_3D_180416']=85
# dataset_ids['CEMBA_5B_180514']=96
# dataset_ids['CEMBA_5B_180529']=97
# dataset_ids['CEMBA_5D_180605']=102
# dataset_ids['CEMBA_5D_180612']=103
# dataset_ids['CEMBA_7B_180423']=115
# dataset_ids['CEMBA_7B_180424']=116

# for dataset in "${!dataset_ids[@]}"
# do
# 	ens_id=${dataset_ids[$dataset]}
# 	# echo $dataset $ens_id
# 	message="Singleton ensemble of $dataset dataset."
# 	ens_name=$dataset
# 	./CEMBA_init_ensemble_v2.py --singleton -d $dataset -ei $ens_id -m $message -en $ens_name
# done

datasets="CEMBA_RS2_Pf10A \
	CEMBA_RS2_Pf10C \
	CEMBA_RS2_Pf11B \
	CEMBA_RS2_Pf12B \
	CEMBA_RS2_Pf3D \
	CEMBA_RS2_Pf4A \
	CEMBA_RS2_Pf5A \
	CEMBA_RS2_Pf6B \
	CEMBA_RS2_Pf7B \
	CEMBA_RS2_Pf9A \
	CEMBA_RS2_Pf9B \
	CEMBA_RS2_Pf9D \
	CEMBA_RS2_Pm10A \
	CEMBA_RS2_Pm10C \
	CEMBA_RS2_Pm11B \
	CEMBA_RS2_Pm12B \
	CEMBA_RS2_Pm3D \
	CEMBA_RS2_Pm4A \
	CEMBA_RS2_Pm5A \
	CEMBA_RS2_Pm6B \
	CEMBA_RS2_Pm7B \
	CEMBA_RS2_Pm9A \
	CEMBA_RS2_Pm9B \
	CEMBA_RS2_Pm9D
"

ens_id=123 # LAST ensemble id!!! 
for dataset in $datasets
do
	ens_id=$((ens_id+1))
	# echo $dataset $ens_id
	message="Singleton ensemble of $dataset dataset."
	ens_name=$dataset
	./CEMBA_init_ensemble_v2.py --singleton -d $dataset -ei $ens_id -m $message -en $ens_name
done

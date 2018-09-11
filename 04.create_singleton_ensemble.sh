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

datasets="MB_v1 \
MB_EA \
MB_EB
"
# 
# 
ens_id=1 # LAST ensemble id!!! 
for dataset in $datasets
do
	ens_id=$((ens_id+1))
	# echo $dataset $ens_id
	message="Singleton ensemble of $dataset dataset."
	ens_name=$dataset
	./CEMBA_init_ensemble_v2.py --singleton -d $dataset -ei $ens_id -m $message -en $ens_name
done

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

# space delimited string
samples="/cndd2/Public_Datasets/CEMBA_cndd2/snmCSeq/Datasets/to_process_upload2.txt"
mapfile -t datasets < $samples
                
ens_id=252  # LAST ensemble id!!! 
for dataset in ${datasets[@]}; do
	ens_id=$((ens_id+1))
	# echo $dataset $ens_id
	message="Singleton ensemble of $dataset dataset."
	ens_name=$dataset
	./CEMBA_init_ensemble_v2.py --singleton -d $dataset -ei $ens_id -m $message -en $ens_name
done

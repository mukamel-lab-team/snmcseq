#!/bin/bash


dataset='CEMBA_SCI_2017'
ens_id=0
message="Singleton ensemble of $dataset dataset."
ens_name=$dataset
./CEMBA_init_ensemble_v2.py --singleton -d $dataset -ei $ens_id -m $message -en $ens_name

# datasets="CEMBA_RS2_Bm3C_rep1 \
# 	CEMBA_RS2_Bm3C_rep2 \
# 	CEMBA_RS2_Bm4B_rep1 \
# 	CEMBA_RS2_Bm4B_rep2 \
# 	CEMBA_RS2_Pf3C \
# 	CEMBA_RS2_Pf4B \
# 	CEMBA_RS2_Pm3C \
# 	CEMBA_RS2_Pm4B \
# 	CEMBA_RS2_Tf3C \
# 	CEMBA_RS2_Tf4B \
# 	CEMBA_RS2_Tf6B \
# 	CEMBA_RS2_Tf7B \
# 	CEMBA_RS2_Tm3C \
# 	CEMBA_RS2_Tm4B \
# 	CEMBA_RS2_Tm6B \
# 	CEMBA_RS2_Tm7B"

# datasets="CEMBA_4C_180417 \
# CEMBA_4C_180419
# "
# 
# 
# ens_id=0 # LAST ensemble id!!! 
# for dataset in $datasets
# do
# 	ens_id=$((ens_id+1))
# 	# echo $dataset $ens_id
# 	message="Singleton ensemble of $dataset dataset."
# 	ens_name=$dataset
# 	./CEMBA_init_ensemble_v2.py --singleton -d $dataset -ei $ens_id -m $message -en $ens_name
# done

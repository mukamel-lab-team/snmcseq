#!/bin/bash

# dataset='CEMBA_3C_171206'
# ens='Ens1'
# message="Singleton ensemble of $dataset dataset."

# dataset='CEMBA_3C_171207'
# ens='Ens2'
# message="Singleton ensemble of $dataset dataset."

# dataset='CEMBA_4B_171212'
# ens='Ens3'
# message="Singleton ensemble of $dataset dataset."

dataset='CEMBA_4B_171213'
ens='Ens4'
message="Singleton ensemble of $dataset dataset."

./CEMBA_init_singleton_ensemble.py -d $dataset -e $ens -m $message
#!/bin/bash

# dataset='CEMBA_3C_171206'
# ens_id=1
# message="Singleton ensemble of $dataset dataset."

# dataset='CEMBA_3C_171207'
# ens_id=2
# message="Singleton ensemble of $dataset dataset."

# dataset='CEMBA_4B_171212'
# ens_id=3
# message="Singleton ensemble of $dataset dataset."

# dataset='CEMBA_4B_171213'
# ens_id=4
# message="Singleton ensemble of $dataset dataset."

dataset='CEMBA_4B_180104'
ens_id=8
message="Singleton ensemble of $dataset dataset."
ens_name="CEMBA_4B_180104"
./CEMBA_init_ensemble_v2.py --singleton -d $dataset -ei $ens_id -m $message -en $ens_name

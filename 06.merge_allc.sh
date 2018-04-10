#!/bin/bash

# ens_id=10
context='CG'
nprocs=2
cluster_type=cluster_mCHmCG_lv_npc50_k5

for ens_id in 10 
do
./CEMBA_merge_allc.py -ei $ens_id -c $context -n $nprocs -cl $cluster_type 
done

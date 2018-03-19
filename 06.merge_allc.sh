#!/bin/bash

# ens_id=7
context='CG'
nprocs=2

for ens_id in 8 9 10 
do
./CEMBA_merge_allc.py -ei $ens_id -c $context -n $nprocs 
done

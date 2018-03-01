#!/bin/bash

# ens_id=7
context='CG'
nprocs=2

for ens_id in 1 2 3 4 5 6
do
./CEMBA_merge_allc.py -ei $ens_id -c $context -n $nprocs 
done

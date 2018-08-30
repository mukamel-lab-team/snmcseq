#!/bin/bash

database='CEMBA'
# datasets="CEMBA_4C_180417 \
# CEMBA_4C_180419
# "
datasets="CEMBA_SCI_2017"

./CEMBA_update_mysql.py -db $database -d $datasets

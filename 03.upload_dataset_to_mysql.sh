#!/bin/bash

database='CEMBA'
datasets="CEMBA_4A_180205 \
CEMBA_4A_180206
"

./CEMBA_update_mysql.py -db $database -d $datasets

#!/bin/bash

database='CEMBA'
datasets="CEMBA_1A_180226 \
CEMBA_1A_180227 \
CEMBA_1C_180208 \
CEMBA_1C_180212
"

./CEMBA_update_mysql.py -db $database -d $datasets

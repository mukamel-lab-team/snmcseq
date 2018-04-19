#!/bin/bash

database='CEMBA'
# datasets="CEMBA_4D_171214 \
# 	CEMBA_4D_171219 \
# 	CEMBA_RS2_Bm3C_rep1 \
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
datasets="CEMBA_3A_180130"

./CEMBA_update_mysql.py -db $database -d $datasets

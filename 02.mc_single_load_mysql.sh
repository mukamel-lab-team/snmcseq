#!/bin/bash
#
# Scans the current directory for files named mc_single_*.tsv and creates corresponding mysql tables
# file format: mc_single_*_$chr.tsv
# table format: mc_single_*_$chr

# Eran, October 2016


db='CEMBA_annoj'

fmysql="mysql -h ocarina -u f7xie -p3405040212"

mkdir -p /scratch/load_mysql_done

$fmysql -e "CREATE DATABASE IF NOT EXISTS $db"
for f in mc_single_*.tsv; do
 echo Loading data from $f into $db
 tab=${f%.tsv};
 $fmysql $db -e "DROP TABLE IF EXISTS $tab"
 $fmysql $db -e "CREATE TABLE IF NOT EXISTS $tab LIKE templates.mc_single"
 $fmysql --local-infile $db -e "LOAD DATA LOCAL INFILE '"$f"' INTO TABLE $tab"
 echo "Done loading $f"
# mv $f /scratch/load_mysql_done
done

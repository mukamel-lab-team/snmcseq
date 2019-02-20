#!/bin/bash
# database, file, table


db=$1
f=$2
tab=$3

fmysql="mysql -h brainome -u f7xie -p3405040212"

# mkdir -p /scratch/load_mysql_done

# $fmysql -e "CREATE DATABASE IF NOT EXISTS $db"
# for f in mc_single_*.tsv; do
#  echo Loading data from $f into $db

 $fmysql $db -e "DROP TABLE IF EXISTS $tab"
 $fmysql $db -e "CREATE TABLE IF NOT EXISTS $tab LIKE templates.mc_single"
 $fmysql --local-infile $db -e "LOAD DATA LOCAL INFILE '"$f"' INTO TABLE $tab"
 
#  echo "Done loading $f"
# mv $f /scratch/load_mysql_done
# done

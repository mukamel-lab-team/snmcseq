#!/bin/bash
prefix=$2
bin=$3
outdir=$4
sort -k2,2n $1 | awk -v bin=$bin -v prefix=$prefix -v outdir=$outdir 'BEGIN {OFS="\t"; assembly=0;pos=0;strand=0;class=0;mc=0;h=0;}
 $1!="chr" {
  if ($2==pos) {
     mc+=$5;h+=$6;
   } else {
     if (assembly!=0) print assembly,pos,strand,class,mc,h > outdir"/allc_"prefix"_"assembly"_"bin ;
    assembly=$1;pos=$2;strand=$3;class=$4;mc=$5;h=$6;
  }
 }
 END { if ($2==pos) print assembly,pos,strand,class,mc,h > outdir"/allc_"prefix"_"assembly"_"bin; }'

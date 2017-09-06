#!/bin/bash
# Function to split an allc table into blocks covering 100kb each
file_in=$1
file_out=`basename $file_in`
mytmpdir=$2
echo "Splitting $file_in and saving as $mytmpdir/$file_out"
zcat -f $file_in | awk -v file_out=$file_out -v mytmpdir=$mytmpdir 'BEGIN {OFS="\t"}
$1!="chr" {
  chr=$1;
  bin=int($2/100000);
  # print $0 > mytmpdir/file_out"_chr"chr"_bin"bin".tsv"
  # print $0 > mytmpdir"/chr"chr"_bin"bin"/"file_out"_chr"chr"_bin"bin".tsv"
  # print $0 > mytmpdir"/chr"chr"_bin"bin"/"file_out
  fn=sprintf("%s/%s_chr%s_bin%0.4d",mytmpdir,file_out,chr,bin)
  print $0 > fn
}'

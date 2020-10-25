#!/bin/bash

pathGrand="/cndd2/Public_Datasets/CEMBA_cndd2/snmCSeq/Datasets/from_ecker_lab/"
pathCurrent=$(pwd)


cd $pathGrand
datasets=($(ls -d [0-9]*))
echo ${datasets[@]}
for dataset in ${datasets[@]}; do
	path="$pathGrand$dataset"
	echo $dataset, $path

	# # rearange allc tables provided by Hanqing
	# path="/cndd2/Public_Datasets/CEMBA_cndd2/snmCSeq/Datasets/from_ecker_lab/10A"
	# basename
	# "allc_CEMBA190423-10A-1-CEMBA190423-10A-2-G8_ad008.tsv.gz"

	cd $path
	for file in $(find . -name "allc_*.gz"); do
		echo $file
		file2="$file.tbi"

		filebase=$(basename $file)
		tmp=${filebase#"allc_CEMBA"}
		tags=(${tmp//"-"/" "})
		date=${tags[0]}
		region=${tags[1]}
		dataset="CEMBA_${region}_${date}"

		fileDst=${filebase/".gz"/".bgz"}
		file2Dst=${filebase/".gz"/".bgz.tbi"}
		mkdir -p "../$dataset/allc"
		mv $file "../$dataset/allc/$fileDst"
		mv $file2 "../$dataset/allc/$file2Dst"
		# break
	done

	# break
done

cd $pathCurrent 

#!/bin/bash

input="/cndd/Public_Datasets/CEMBA/snmCSeq/Datasets/CEMBA_RS2_Pf10A/allc \
	/cndd/Public_Datasets/CEMBA/snmCSeq/Datasets/CEMBA_RS2_Pf10C/allc \
	/cndd/Public_Datasets/CEMBA/snmCSeq/Datasets/CEMBA_RS2_Pf11B/allc \
	/cndd/Public_Datasets/CEMBA/snmCSeq/Datasets/CEMBA_RS2_Pf12B/allc \
	/cndd/Public_Datasets/CEMBA/snmCSeq/Datasets/CEMBA_RS2_Pf3D/allc \
	/cndd/Public_Datasets/CEMBA/snmCSeq/Datasets/CEMBA_RS2_Pf4A/allc \
	/cndd/Public_Datasets/CEMBA/snmCSeq/Datasets/CEMBA_RS2_Pf5A/allc \
	/cndd/Public_Datasets/CEMBA/snmCSeq/Datasets/CEMBA_RS2_Pf6B/allc \
	/cndd/Public_Datasets/CEMBA/snmCSeq/Datasets/CEMBA_RS2_Pf7B/allc \
	/cndd/Public_Datasets/CEMBA/snmCSeq/Datasets/CEMBA_RS2_Pf9A/allc \
	/cndd/Public_Datasets/CEMBA/snmCSeq/Datasets/CEMBA_RS2_Pf9B/allc \
	/cndd/Public_Datasets/CEMBA/snmCSeq/Datasets/CEMBA_RS2_Pf9D/allc \
	/cndd/Public_Datasets/CEMBA/snmCSeq/Datasets/CEMBA_RS2_Pm10A/allc \
	/cndd/Public_Datasets/CEMBA/snmCSeq/Datasets/CEMBA_RS2_Pm10C/allc \
	/cndd/Public_Datasets/CEMBA/snmCSeq/Datasets/CEMBA_RS2_Pm11B/allc \
	/cndd/Public_Datasets/CEMBA/snmCSeq/Datasets/CEMBA_RS2_Pm12B/allc \
	/cndd/Public_Datasets/CEMBA/snmCSeq/Datasets/CEMBA_RS2_Pm3D/allc \
	/cndd/Public_Datasets/CEMBA/snmCSeq/Datasets/CEMBA_RS2_Pm4A/allc \
	/cndd/Public_Datasets/CEMBA/snmCSeq/Datasets/CEMBA_RS2_Pm5A/allc \
	/cndd/Public_Datasets/CEMBA/snmCSeq/Datasets/CEMBA_RS2_Pm6B/allc \
	/cndd/Public_Datasets/CEMBA/snmCSeq/Datasets/CEMBA_RS2_Pm7B/allc \
	/cndd/Public_Datasets/CEMBA/snmCSeq/Datasets/CEMBA_RS2_Pm9A/allc \
	/cndd/Public_Datasets/CEMBA/snmCSeq/Datasets/CEMBA_RS2_Pm9B/allc \
	/cndd/Public_Datasets/CEMBA/snmCSeq/Datasets/CEMBA_RS2_Pm9D/allc
"


./CEMBA_run_mc_gene_level.py -f -i $input -n 2 


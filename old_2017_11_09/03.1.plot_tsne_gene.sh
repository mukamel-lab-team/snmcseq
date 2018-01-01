#!/bin/bash

gene_dir=/cndd/Public_Datasets/single_cell_methylome/gene_level/human/genebody_MB_EB_by_gene

for gene in GAD1 SATB2 TYRO3 ARPP21 SLC17A7 TBR1 CAMK2A ITPKA CUX1 CUX2 RORB DEPTOR VAT1L SULF1 TLE4 FOXP2 GRIK3 BCL6 ERBB4 GAD1 SLC6A1 ADARB2 PROX1 SV2C PVALB SOX6 RELN CACNA2D2 LHX6 GRIA1
# for gene in GAD1
do
#	grep -i $gene ${gene_dir}/Pool_2256_AD010_indexed_R1_mch_genebody.txt | wc -l
 	for perplexity in 25
 	do
 		input=./tsne/human_combined_MB_EB_CH_100000_out_v1_${perplexity}.tsv
 		output=./tsne/marker_genes/$(basename $input | cut -f 1 -d '.')_${gene}_normalized.pdf
 
 		./scripts_plot/plot_tsne_v2.py -i $input --gene $gene --gene_dir $gene_dir -o $output --title $gene'_perplexity='${perplexity}'_normalized_mCH' --normalize 
 	done 
done

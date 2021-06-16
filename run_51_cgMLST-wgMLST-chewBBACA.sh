#!/bin/bash
###############################################################################
#51 cgMLST and wgMLST

###############################################################################

echo "started.... step-51 cgMLST and wgMLST ----------------------------------"

###############################################################################

if [ ! -d results/51_cg-wg-MLST ] ; then
mkdir results/51_cg-wg-MLST
fi


## step-0 Create training file

prodigal -i closed_genomes_fna/NC_003210.fna -t Listeria_monocytogenes.trn -p single

## step -1 Create Scheme

chewBBACA.py CreateSchema -i closed_genomes/ --cpu 6 -o schema_seed --ptf Streptococcus_agalactiae.trn

## step -2 AlleleCall

chewBBACA.py AlleleCall -i test_genomes/ -g schema_seed/ -o AlleleCall_results --cpu 6 -b /home/swapnil/tools/ncbi-blast-2.7.1+/bin/blastp --ptf Streptococcus_agalactiae.trn

## step -3 TestGenomeQuality

chewBBACA.py TestGenomeQuality -i results_cg2/results_20180426T150016/alleles.tsv -n 12 -t 200 -s 5 -o OutFolder

## step -4 ExtractCgMLST

chewBBACA.py ExtractCgMLST -i results_cg2/results_20180426T150016/results_alleles.tsv -o output_folders -p 0.99

## step -5 SchemaEvaluator

chewBBACA.py SchemaEvaluator -i schema_seed/ -ta 11 -l test.html --cpu 6 --title "test"

#open html by firefox
















































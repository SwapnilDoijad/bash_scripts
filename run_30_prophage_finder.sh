#!/bin/bash
###############################################################################
#30 Prophage finder

###############################################################################

echo "started.... prophage finder --------------------------------------------"

###############################################################################

if [ ! -d results/30_prophage_finder ]; then
mkdir results/30_prophage_finder
fi

if [ ! -d results/30_prophage_finder/gff3 ]; then
mkdir results/30_prophage_finder/gff3
fi

###############################################################################

for F1 in $(cat list.txt);do

echo "running prophage finder for $F1"

(bp_genbank2gff3 results/08_annotation/all_gbk/$F1.gbk -o results/30_prophage_finder/gff3/) > /dev/null 2>&1
sed -i '/FASTA/,$d' results/30_prophage_finder/gff3/$F1.gbk.gff
sed -i '/ssrA/d' results/30_prophage_finder/gff3/$F1.gbk.gff
(perl /home/swapnil/pipeline/tools/ProphET-master/ProphET_standalone.pl --fasta_in results/04_assembly/all_fasta/$F1.fasta --gff_in results/30_prophage_finder/gff3/$F1.gbk.gff --outdir results/30_prophage_finder/$F1) > /dev/null 2>&1

echo "finished prophage finder for $F1"

done


###############################################################################

echo "Completed.... prophage finder ------------------------------------------"

###############################################################################


#!/bin/bash
###############################################################################
#24A plasmid finder Plasflow

###############################################################################

echo "started.... Plasmid search ---------------------------------------------"

###############################################################################

echo "provide list file (for e.g. all)"
read l
list=$(echo "list.$l.txt")

if [ ! -d results/24A_plasFlow_"$l" ]; then
mkdir results/24A_plasFlow_"$l"
fi

for F1 in $(cat $list); do
mkdir results/24A_plasFlow_"$l"/$F1
python /home/swapnil/tools/PlasFlow/PlasFlow.py --input results/04_assembly/raw_files/$F1/contigs.fasta --output $F1
mv *.fasta results/24A_plasFlow_"$l"/$F1/
mv $F1 results/24A_plasFlow_"$l"/$F1/
done


###############################################################################

echo "Completed.. Plasmid search ---------------------------------------------"

###############################################################################

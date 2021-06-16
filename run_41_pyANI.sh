#!/bin/bash
###############################################################################
#41 pyANI
# files should end with fasta(?)

###############################################################################

echo "started.... step-41 pyANI ----------------------------------------------"

###############################################################################
if [ ! -d results/41_pyani ] ; then
mkdir results/41_pyani
fi

if [ ! -d results/41_pyani/fasta ] ; then
mkdir results/41_pyani/fasta
fi

#------------------

for F1 in $(cat list.txt); do
cp results/04_assembly/all_fasta/$F1.fasta results/41_pyani/fasta
done

#------------------

average_nucleotide_identity.py -i results/41_pyani/fasta -o results/41_pyani/results -m ANIblastall -g

rm results/41_pyani/tmp

#------------------
# get a representative of a clutster at 95% (-h 5)

echo "running ggrasp.R"

/home/swapnil/tools/GGRaSP-master/ggrasp.R -i results/41_pyani/results/ANIblastall_percentage_identity.tab -d 100 -h 5 -o results/41_pyani/representative_of_a_clusters_at_95_ANI --plottree 

rm Rplots.pdf

###############################################################################

echo "completed.... step-41 pyANI --------------------------------------------"

###############################################################################

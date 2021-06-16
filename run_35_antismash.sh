#!/bin/bash
###############################################################################
## check .bashrc file and remove # , see below

#------------------
# added by Miniconda2 installer
# export PATH="/home/swapnil/miniconda2/bin:$PATH"

# added by Anaconda2 installer
# export PATH="/home/swapnil/anaconda2/bin:$PATH"
#------------------

###############################################################################
echo "started.... antiSMASH ----------------------------------------------------"
###############################################################################
    echo "provide list file (for e.g. all)"
    echo "-------------------------------------------------------------------------"
    ls list.*.txt | sed 's/ /\n/g'
    echo "-------------------------------------------------------------------------"
    read l
    list=$(echo "list.$l.txt")

(mkdir results/35_antiSMASH) > /dev/null 2>&1 
(mkdir results/35_antiSMASH/results) > /dev/null 2>&1 
###############################################################################

source activate antismash

for F1 in $(cat $list); do

    echo "running antiSMASH $F1"

    antismash results/08_annotation/all_gbk/$F1.gbk

    mv results/08_annotation/all_gbk/$F1 results/35_antiSMASH/results/

    echo "finished antiSMASH $F1"

done

source deactivate antismash

###############################################################################
echo "completed.... antiSMASH --------------------------------------------------"
###############################################################################

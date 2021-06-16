#!/bin/bash

## in conda my_env_3.6.3 environment
###############################################################################
# MLST
###############################################################################
echo "started...... step-15 MLST ---------------------------------------------"
###############################################################################
## initial file and directory preparation
    source /home/swapnil/miniconda3/etc/profile.d/conda.sh
    conda activate my_env_3.6.3

    echo "provide list file (for e.g. all)"
    echo "-------------------------------------------------------------------------"
    ls list.*.txt | sed 's/ /\n/g'
    echo "-------------------------------------------------------------------------"
    read l
    list=$(echo "list.$l.txt")

    echo "provide species (without space for e.g. cloacae)"
    read species

    echo "choose the scheme"
    grep $species /home/swapnil/tools/mlst-master/db/species_scheme_map.tab | awk '{print $NF}'

    read scheme 

    (mkdir results/15_MLST) > /dev/null 2>&1
    (mkdir results/15_MLST/fasta) > /dev/null 2>&1
###############################################################################
## run MLST
    for F1 in $(cat $list) ; do
    cp results/04_assembly/all_fasta/$F1.fasta results/15_MLST/fasta/$F1.fasta
    done

    echo "running MLST"
    mlst --quiet --threads=8 --nopath --scheme="$scheme" --novel=results/15_MLST/novel-allels.fa results/15_MLST/fasta/*.fasta >> results/15_MLST/mlst.tab
    ## without species description / automatic mlst
    # /home/swapnil/pipeline/tools/mlst-master/bin/mlst --quiet --threads=8 --nopath --novel=results/15_MLST/novel-allels.fa results/04_assembly/all_fasta/*.fasta >> results/15_MLST/mlst.tab
    rm -rf results/15_MLST/fasta
    conda deactivate
###############################################################################
echo "started...... step-15 MLST ---------------------------------------------"
###############################################################################


#!/bin/bash
###############################################################################
    source /home/swapnil/miniconda3/etc/profile.d/conda.sh
    conda activate myenv
###############################################################################

    echo "provide list file (for e.g. all)"
    echo "---------------------------------------------------------------------"
    ls list.*.txt | sed 's/ /\n/g' | sed 's/list\.//g' | sed 's/\.txt//g'
    echo "---------------------------------------------------------------------"
    read l
    list=$(echo "list.$l.txt")

    (mkdir results) > /dev/null 2>&1
    (mkdir results/06_identification_gtdb) > /dev/null 2>&1
    (mkdir results/06_identification_gtdb/fasta) > /dev/null 2>&1
    (mkdir results/06_identification_gtdb/results) > /dev/null 2>&1

###############################################################################
## run gtdbtk
    ## copy files in a directory
    for F1 in $(cat $list);do
        echo "copying $F1"
        cp results/04_assembly/all_fasta/$F1.fasta results/06_identification_gtdb/fasta/
    done
    
    ## running gtdb-tk
    echo "running gtdb-tk"
    gtdbtk classify_wf --cpus 12 --genome_dir results/06_identification_gtdb/fasta/ --out_dir results/06_identification_gtdb/results --extension fasta

    ## rm -rf results/06_identification_gtdb/fasta
###############################################################################
conda deactivate myenv
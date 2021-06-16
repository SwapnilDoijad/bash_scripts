#!/bin/bash
###############################################################################
# 12 harvest suite
###############################################################################
echo "started...... step-11 harvest suite ------------------------------------"
###############################################################################
    echo "provide list file (for e.g. all)"
    echo "---------------------------------------------------------------------"
    ls list.*.txt | sed 's/ /\n/g' | sed 's/list\.//g' | sed 's/\.txt//g'
    echo "---------------------------------------------------------------------"
    read l
    list=$(echo "list.$l.txt")

    (mkdir results/12_parsnp_"$l") > /dev/null 2>&1
    (mkdir results/12_parsnp_"$l"/fasta) > /dev/null 2>&1

###############################################################################
## copy files
    echo "copying files"
    for F1 in $(cat $list); do
        cp results/04_assembly/all_fasta/$F1.fasta results/12_parsnp_"$l"/fasta/
    done

## run parsnp
echo "running parsnp program"
    parsnp -g results/00_ref/ref.*.gbk -c -d results/12_parsnp_"$l"/fasta/ -p 16 -o results/12_parsnp_"$l"/ 

## 
    echo "parnsp finished, cleaning"
    (rm -rf results/12_parsnp_"$l"/fasta/) > /dev/null 2>&1
###############################################################################
echo "completed.... step-11 harvest suite ------------------------------------"
###############################################################################

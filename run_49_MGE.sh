#!/bin/bash
###############################################################################
# 49 MEF
###############################################################################
echo "started... step-49 Mobile element finder -------------------------------------------"
###############################################################################
    echo "provide list file (for e.g. all)"
    echo "---------------------------------------------------------------------"
    ls list.*.txt | awk -F'.' '{print $2}'
    echo "---------------------------------------------------------------------"
    read l
    list=$(echo "list.$l.txt")

    (mkdir results/49_MEF) > /dev/null 2>&1
    (mkdir results/49_MEF/raw_files) > /dev/null 2>&1
###############################################################################
for F1 in $(cat $list); do

    echo "running.... step-49 Mobile element finder for..." $F1
    (mefinder find --contig results/04_assembly/all_fasta/"$F1".fasta results/49_MEF/raw_files/$F1) > /dev/null 2>&1
    echo "finished... step-49 Mobile element finder for..." $F1

done


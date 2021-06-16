#!/bin/bash
###############################################################################
# 43 Island finder
###############################################################################
echo "started... step-43 Island finder ---------------------------------------"
###############################################################################
## intial file and directory preparations
    echo "provide list file (for e.g. all)"
    echo "---------------------------------------------------------------------"
    ls list.*.txt | sed 's/ /\n/g'
    echo "---------------------------------------------------------------------"
    read l
    list=$(echo "list.$l.txt")
    
    read -p "please be out of conda environment, if any"

(mkdir results/43_island_finder) > /dev/null 2>&1
(mkdir results/43_island_finder/raw_result_files) > /dev/null 2>&1
###############################################################################
for F1 in $(cat $list); do
    if [ -d results/43_island_finder/raw_result_files/$F1 ] ; then
    echo "$F1 already finished"
    else
        echo "running.... step-43 island finder... $F1"
        (mkdir results/43_island_finder/raw_result_files/$F1) > /dev/null 2>&1
        #(
        perl /home/swapnil/tools/Islander_software/Islander.pl results/08_annotation/raw_files/$F1/$F1 --verbose --translate --trna --annotate --reisland --table 11 --nocheck #) > /dev/null 2>&1
        mv all_final_arrayed_islands.gff results/43_island_finder/raw_result_files/$F1/$F1.islands.gff
        mv check_tandem.txt results/43_island_finder/raw_result_files/$F1/$F1.check_tandem.txt
    fi
done
###############################################################################
echo "finished... step-43 Island finder --------------------------------------"
###############################################################################

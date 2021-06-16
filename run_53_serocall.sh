#!/bin/bash
###############################################################################
#53 SeroCall (for SeroCall can identify and quantitate the capsular serotypes 
#   in Illumina whole-genome sequencing samples of S. pneumoniae)
###############################################################################
    echo "provide list file (for e.g. all)"
    echo "-------------------------------------------------------------------------"
    ls list.*.txt | sed 's/ /\n/g'
    echo "-------------------------------------------------------------------------"
    read l
    list=$(echo "list.$l.txt")

    echo "raw data path? (for e.g. /media/network/reads_database/p_Entb_Germany)"
    read path

    (mkdir $path/results) > /dev/null 2>&1
    (mkdir $path/results/53_serocall) > /dev/null 2>&1
    (mkdir $path/results/53_serocall/raw_files) > /dev/null 2>&1
###############################################################################

for F1 in $(cat $list); do
    (mkdir $path/results/53_serocall/raw_files/$F1) > /dev/null 2>&1
    /home/swapnil/tools/SeroCall-master/serocall \
    $path/results/02_filtered_reads/$F1/"$F1"_R1_001.filtered_paired.fastq.gz \
    $path/results/02_filtered_reads/$F1/"$F1"_R2_001.filtered_paired.fastq.gz \
    -o $path/results/53_serocall/raw_files/$F1/$F1
    -t 8
done 
###############################################################################
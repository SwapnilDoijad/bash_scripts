#!/bin/bash
###############################################################################
#18 scoary
#NOTE: if the strain name contains only numeric characters, there is some problem running scoary
#NOTE: manually remove comma , from gene_presence_absence.csv file
#NOte: prepare traits.csv file, opne with LibreOffice Calc and again save (in csv format)
###############################################################################

echo "started.... scoary -----------------------------------------------------"

###############################################################################
## initial inout and file preparation
    echo "Which folder to look for gene_presence_absence.csv?"
    echo "---------------------------------------------------------------------"
    ls list.*.txt | sed 's/ /\n/g' | sed 's/list\.//g' | sed 's/\.txt//g'
    echo "---------------------------------------------------------------------"
    read s1
    s2=$(echo "_$s1")

    echo "provide trait file (for e.g. all)"
    echo "---------------------------------------------------------------------"
    ls data/traits.scoary.*.csv | sed 's/data\/traits\.scoary\.//g' | sed 's/\.csv//g'
    echo "---------------------------------------------------------------------"
    read traits_tmp
    traits=$(echo "data/traits.scoary.$traits_tmp.csv" )

    echo "want to suffix files? type "y" else press enter to continue"
    read answer
    if [ "$answer" = "y" ]; then
    s=$(echo "_$traits_tmp")
    fi

    mkdir results/18_scoary"$s"
    mkdir results/18_scoary"$s"/tmp

    cp results/17_roary"$s2"/roary_results/gene_presence_absence.csv results/18_scoary"$s"/tmp/

    sed -i 's/ /_/g' results/18_scoary"$s"/tmp/gene_presence_absence.csv
   # sed -i 's/,/_/g' results/18_scoary"$s"/tmp/gene_presence_absence.csv

###############################################################################

scoary -g results/18_scoary"$s"/tmp/gene_presence_absence.csv -t $traits -o results/18_scoary"$s"/

ssconvert results/18_scoary"$s"/*.csv results/18_scoary"$s"/results.xls
mv results/18_scoary"$s"/*.csv results/18_scoary"$s"/tmp/

rm results/18_scoary"$s"/*.log

###############################################################################

echo "completed.... scoary ---------------------------------------------------"

###############################################################################

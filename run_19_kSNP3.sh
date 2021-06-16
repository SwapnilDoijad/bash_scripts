#!/bin/bash
###############################################################################
#19 kSNP3

###############################################################################

echo "started.... kSNP3 ------------------------------------------------------"

###############################################################################

    echo "provide list file (for e.g. all)"
    echo "---------------------------------------------------------------------"
    ls list.*.txt | sed 's/ /\n/g'
    echo "---------------------------------------------------------------------"
    read l
    list=$(echo "list.$l.txt")

    s=$(echo "_$l")
    
    (mkdir results/19_kSNP3"$s") > /dev/null 2>&1
    (mkdir results/19_kSNP3"$s"/tmp) > /dev/null 2>&1
    (mkdir results/19_kSNP3"$s"/fasta) > /dev/null 2>&1
###############################################################################

    J=1
    for F1 in $(cat $list); do
    cp results/04_assembly/all_fasta/"$F1".fasta results/19_kSNP3"$s"/fasta/strain"$J".fasta; 
    echo strain"$J" $F1 >> results/19_kSNP3"$s"/tmp/re-labeling.tab; 
    let J=J+1
    done

    sed -i 's/ /\t/g' results/19_kSNP3"$s"/tmp/re-labeling.tab ;
    sed -i '1 i\#old_name	new_name' results/19_kSNP3"$s"/tmp/re-labeling.tab

    cd results/19_kSNP3"$s"/fasta/
    ls *.fasta > fasta_list.txt
    sed -i 's/\.fasta//g' fasta_list.txt
    cd ..
    cd ..
    cd ..

#------------------------------------------------------------------------------
#KChooser

    echo "running Kchooser"

    cat results/19_kSNP3"$s"/fasta/*.fasta > results/19_kSNP3"$s"/tmp/all.fasta

    (Kchooser results/19_kSNP3"$s"/tmp/all.fasta) > /dev/null 2>&1

    grep "The optimum value of K is " Kchooser.report | sed 's/The optimum value of K is//g' | sed 's/\.//g' | sed 's/ //g' > results/19_kSNP3"$s"/tmp/kmer.txt

    rm jellyout.txt
    rm Kchooser.report
    rm medSeq.fasta
    rm shortSeq.fasta
    rm results/19_kSNP3"$s"/tmp/all.fasta

    echo "finished Kchooser"

    #------------------------------------------------------------------------------
    (rm results/19_kSNP3"$s"/tmp/in_file.txt)> /dev/null 2>&1
    WorkDir=$(echo $PWD)
    for F2 in $(cat results/19_kSNP3"$s"/fasta/fasta_list.txt); do
    echo "$WorkDir/results/19_kSNP3"$s"/fasta/$F2.fasta" $F2 >> results/19_kSNP3"$s"/tmp/in_file.txt	######### /home/swapnil/pipeline/
    sed -i 's/ /\t/g' results/19_kSNP3"$s"/tmp/in_file.txt
    done

    V2=$(cat results/19_kSNP3"$s"/tmp/kmer.txt)

    kSNP3 -in results/19_kSNP3"$s"/tmp/in_file.txt -k 11 -outdir results/19_kSNP3"$s"/result

    #kSNP3 -in results/19_kSNP3"$s"/tmp/in_file.txt -k $V2 -outdir results/19_kSNP3"$s"/results -core -min_frac 0.75

###############################################################################

echo "completed.... kSNP3 ----------------------------------------------------"

###############################################################################

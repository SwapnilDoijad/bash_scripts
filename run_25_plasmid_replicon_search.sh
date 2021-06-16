#!/bin/bash
###############################################################################
#25 plasmid serach complete
###############################################################################
echo "started.... plasmid search complete -------------------------------------"
###############################################################################
## file and directory preparation
    echo "provide list file (for e.g. all)"
    echo "---------------------------------------------------------------------"
    ls list.*.txt | sed 's/ /\n/g'
    echo "---------------------------------------------------------------------"
    read l
    list=$(echo "list.$l.txt")
    s=$(echo "_$l")

    (mkdir results/25_plasmid_replicons"$s")> /dev/null 2>&1
    (mkdir results/25_plasmid_replicons"$s"/tmp)> /dev/null 2>&1
    (mkdir results/25_plasmid_replicons"$s"/raw_files)> /dev/null 2>&1
###############################################################################

(rm results/25_plasmid_replicons"$s"/raw_files/$F1/tmp/$F1.$WGS_replicon.list)> /dev/null 2>&1
for F1 in $(cat $list); do
    echo "runnning... Plasmid search for $F1---------------------------------------"
    blastn -db /media/swapnil/share/databases/plasmid/plasmid.fsa -query results/04_assembly/all_fasta/$F1.fasta -max_hsps 3 -evalue 1e-100 -num_threads 8 -outfmt "6 sseqid qseqid sstart send qstart qend slen qlen evalue bitscore length mismatch gaps pident qcovs" > results/25_plasmid_replicons"$s"/raw_files/$F1.plasmid.blast.tmp
    WGS_pReplicon=$(cat results/25_plasmid_replicons"$s"/raw_files/$F1.plasmid.blast.tmp | awk -F'_' '{print $1}' | sort -u | tr "\n" "_" | sed '$s/.$//')
    if [ -z "$WGS_pReplicon" ]; then 
        WGS_pReplicon=$(echo "NA")
        else
        cat results/25_plasmid_replicons"$s"/raw_files/$F1.plasmid.blast.tmp | awk '{print $2}' | sed -e "s/$/\.fa/" > results/25_plasmid_replicons"$s"/raw_files/$F1.WGS_pReplicon.contigs.list
    fi
    echo $F1 $WGS_pReplicon > results/25_plasmid_replicons"$s"/raw_files/$F1.WGS_pReplicon.list
    echo $F1 $WGS_pReplicon >> results/25_plasmid_replicons"$s"/All.WGS_pReplicon.csv
done
###############################################################################
echo "Completed.. plasmid serach complete ------------------------------------"
###############################################################################

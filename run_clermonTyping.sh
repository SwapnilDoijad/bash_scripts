
###############################################################################
    echo "provide list file (for e.g. all)"
    echo "---------------------------------------------------------------------"
    ls list.*.txt | awk -F'.' '{print $2}'
    echo "---------------------------------------------------------------------"
    read l
    list=$(echo "list.$l.txt")

    (mkdir results/clermonTyping) > /dev/null 2>&1
    (mkdir results/clermonTyping/raw_files) > /dev/null 2>&1
###############################################################################

for F1 in $(cat $list); do
echo "ClermonTyping running for $F1"
(/home/swapnil/tools/ClermonTyping-master/clermonTyping.sh --fasta results/04_assembly/all_fasta/$F1.fasta) > /dev/null 2>&1
mv analysis*/ results/clermonTyping/raw_files/$F1

path1=$(find results/clermonTyping/raw_files/$F1/ -name '*phylogroups.txt')
CT=$(cat $path1 | awk -F'\t' '{print $5}' )
echo $F1 $CT >> results/clermonTyping/results.csv
done

###############################################################################

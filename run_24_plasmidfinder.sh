###############################################################################
# 08 annotation
###############################################################################
echo "started... step-08 annotation -------------------------------------------"
###############################################################################
    echo "provide list file (for e.g. all)"
    echo "---------------------------------------------------------------------"
    ls list.*.txt | awk -F'.' '{print $2}'
    echo "---------------------------------------------------------------------"
    read l
    list=$(echo "list.$l.txt")

    (mkdir results/plasmidfinder) > /dev/null 2>&1
    (mkdir results/plasmidfinder/raw_files) > /dev/null 2>&1
###############################################################################
(rm results/plasmidfinder/results.csv) > /dev/null 2>&1

for F1 in $(cat $list); do
##contigs
    echo "PlasmidFinder running for $F1"
    (mkdir results/plasmidfinder/raw_files/$F1) > /dev/null 2>&1
    (/home/swapnil/tools/plasmidfinder/plasmidfinder.py \
    -i results/04_assembly/"Escherichia_coli_"$F1.fa \
    -o results/plasmidfinder \
    -tmp results/plasmidfinder/raw_files/$F1 \
    -p /media/swapnil/databases/bioinfoDBs/plasmidfinder_db) > /dev/null 2>&1

    mv results/plasmidfinder/data.json results/plasmidfinder/raw_files/$F1/$F1.data.json

    V1=$(jq . results/plasmidfinder/raw_files/$F1/$F1.data.json | grep 'plasmid' | sed '1d' | sed 's/"//g' | sed 's/,//g' | sed 's/plasmid://g' | tr '\r\n' ' ' | sed 's/             //g' | sed 's/ /,/g' | sed 's/.\{1\}$/\n/')
    echo "$F1 $V1" >> results/plasmidfinder/results.csv
done

exit


## fastq
    /home/swapnil/tools/plasmidfinder/plasmidfinder.py \
    -i /media/swapnil/network/reads_database/p_Enterobacter_Germany/results/02_filtered_reads/RBF-18-0058-1/RBF-18-0058-1_R1_001.filtered_paired.fastq.gz \
    -o plasmid_results \
    -tmp tmp_dir\
    -mp /home/swapnil/tools/kma/kma \
    -p /media/swapnil/databases/bioinfoDBs/plasmidfinder_db
## parse json results

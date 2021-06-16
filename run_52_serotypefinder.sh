###############################################################################
    echo "provide list file (for e.g. all)"
    echo "---------------------------------------------------------------------"
    ls list.*.txt | awk -F'.' '{print $2}'
    echo "---------------------------------------------------------------------"
    read l
    list=$(echo "list.$l.txt")

    (mkdir results/52_SerotypeFinder) > /dev/null 2>&1
    (mkdir results/52_SerotypeFinder/raw_files) > /dev/null 2>&1
###############################################################################
STFinder_DB=/media/swapnil/databases/bioinfoDBs/serotypefinder_db
(rm results/52_SerotypeFinder/results.csv) > /dev/null 2>&1
for F1 in $(cat $list); do
echo "SerotypeFinder running for $F1"
#(mkdir results/52_SerotypeFinder/raw_files/$F1) > /dev/null 2>&1
#(docker run --rm -it -v $STFinder_DB:/database -v $(pwd):/workdir serotypefinder -i results/04_assembly/all_fasta/$F1.fasta -o results/52_SerotypeFinder/raw_files/$F1) > /dev/null 2>&1
V1=$(cat results/52_SerotypeFinder/raw_files/$F1/data.json | tr , '\n' | tr { '\n' | grep serotype | awk -F'"' 'FNR>1 {print $4}' | tr '\n' ' ' )
echo $V1
O_type=$(echo $V1 | tr ' ' '\n' | grep O | head -1)
H_type=$(echo $V1 | tr ' ' '\n' | grep H | head -1)
echo $O_type 
echo $H_type
echo $F1 $O_type $H_type >> results/52_SerotypeFinder/results.csv
done
###############################################################################

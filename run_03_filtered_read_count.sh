#!/bin/bash
###############################################################################
#3 filtered read count

###############################################################################

echo "started...... step-3 filtered read count -------------------------------"

###############################################################################
echo "provide list file (for e.g. all)"
ls list.*.txt | sed 's/ /\n/g'
read l
list=$(echo "list.$l.txt")

echo "raw data path? (for e.g. /media/network/reads_database/p_Entb_Germany)"
ls /media/network/reads_database/
read path_tmp1
path=$(echo "/media/network/reads_database/$path_tmp1")

(mkdir results/03_filtered_reads_count) > /dev/null 2>&1 
##-----------------------------------------------------------------------------

for F1 in $(cat $list); do

if [ ! -f $path/results/03_filtered_reads_count/$F1/"$F1".3_filtered_read_count.statistics.tab ]; then
##-------
## waiting for the file
while [ ! -f $path/results/02_filtered_reads/$F1/"$F1"_R1_001.filtered_paired.fastq.gz ] ; do
sleep 1m 
done
##-------

echo "running.... step-3 filtered_read_count for..." $F1

if [ ! -d results/03_filtered_reads_count/$F1 ]; then
mkdir results/03_filtered_reads_count/$F1
fi

if [ ! -d results/03_filtered_reads_count/$F1/tmp ]; then
mkdir results/03_filtered_reads_count/$F1/tmp
fi

V1=$(zcat $path/results/02_filtered_reads/$F1/"$F1"_R1_001.filtered_paired.fastq.gz | echo $((`wc -l`/4)))
V2=$(zcat $path/results/02_filtered_reads/$F1/"$F1"_R2_001.filtered_paired.fastq.gz | echo $((`wc -l`/4)))

zcat $path/results/02_filtered_reads/$F1/"$F1"_R1_001.filtered_paired.fastq.gz | awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' > results/03_filtered_reads_count/$F1/tmp/$F1.3_filtered_read_count.tmp.statistics.tab
V3=$(cat results/03_filtered_reads_count/$F1/tmp/$F1.3_filtered_read_count.tmp.statistics.tab | awk '{$3=$1*$2; C+=$3; B+=$2; ARL= C / B} END {print ARL}') 
rm results/03_filtered_reads_count/$F1/tmp/$F1.3_filtered_read_count.tmp.statistics.tab

zcat $path/results/02_filtered_reads/$F1/"$F1"_R2_001.filtered_paired.fastq.gz | awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' > results/03_filtered_reads_count/$F1/tmp/$F1.3_filtered_read_count.tmp.statistics.tab
V4=$(cat results/03_filtered_reads_count/$F1/tmp/$F1.3_filtered_read_count.tmp.statistics.tab | awk '{$3=$1*$2; C+=$3; B+=$2; ARL= C / B} END {print ARL}') 
rm results/03_filtered_reads_count/$F1/tmp/$F1.3_filtered_read_count.tmp.statistics.tab

V5=$(awk "BEGIN {print ($V3+$V4)/2; exit}")

echo $F1 $V1 $V2 $V5 > results/03_filtered_reads_count/$F1/tmp/$F1.3_filtered_read_count.statistics.tab

(V6=$(grep "$F1" results/03_filtered_reads_count/filtered_read_count.statistics.tab)) > /dev/null 2>&1 
if [ "$V6" !=  "$F1" ]; then
echo $F1 $V1 $V2 $V5 >> results/03_filtered_reads_count/filtered_read_count.statistics.tab
fi

if [ ! -f results/03_filtered_reads_count/$F1/tmp/$F1.3_filtered_read_count.statistics.tab ]; then
echo "failed! $F1 filtered read count " 
echo "failed! $F1 filtered read count " >> results/failed.txt
fi

echo "finished... step-3 filtered_read_count for..." $F1
else
echo "Reads already filtered for $F1 "
fi
done

###############################################################################

echo "completed.. step-3 filtered read count ---------------------------------"

###############################################################################

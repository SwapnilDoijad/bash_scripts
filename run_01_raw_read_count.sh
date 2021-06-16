#!/bin/bash
###############################################################################
#1 raw read count
## Warning: Do not save orignial four files of NextSeq  and a single merged file *R1*.fastq.gz in one folder, else the read count will be performed for both
#NanoStat --fastq /media/network/reads_database/p_Entb_Germany/"$ncbi"data/nanopore/final_reads/Survcare220_BC12.fastq.gz --outdir /media/swapnil/Swapnil_Sony_HD-E1/"$ncbi"data/Genome_database/p_Eb_two_Germany/results/01_raw_read_count/00_nanopore -p Survcare220_BC12 -n Survcare220_BC12
###############################################################################
echo "started...... step-1 raw_read_count ------------------------------------"
###############################################################################
echo "provide list file (for e.g. all)"
ls list.*.txt | sed 's/ /\n/g'
read l
list=$(echo "list.$l.txt")

echo "raw data path? (for e.g. /media/network/reads_database/p_Entb_Germany)"
ls /media/network/reads_database/
read path_tmp1
path=$(echo "/media/network/reads_database/$path_tmp1")

echo "do you want to analyze the NCBI data? type y for yes" 
read ncbi_tmp
ncbi=$(echo "ncbi/")

(mkdir $path/"$ncbi"data) > /dev/null 2>&1 
(mkdir $path/"$ncbi"data/illumina) > /dev/null 2>&1 
(mkdir $path/"$ncbi"data/illumina/final_raw_reads) > /dev/null 2>&1 

(rm results/01_raw_read_count/raw_read_count.statistics.tab) > /dev/null 2>&1 
for F1 in $(cat $list);do

##-----------------------------------------------------------------------------
## Identify the raw reads of Miseq/NextSeq and place them in final_raw_reads

V1=$(ls $path/"$ncbi"data/illumina/raw_reads/"$F1"*.gz | wc -l )
if (( $V1 > 2 )); then
(cat $path/"$ncbi"data/illumina/raw_reads/"$F1"_*R1*.gz > $path/"$ncbi"data/illumina/final_raw_reads/"$F1"_merged_R1.fastq.gz) > /dev/null 2>&1 
(cat $path/"$ncbi"data/illumina/raw_reads/"$F1"_*R2*.gz > $path/"$ncbi"data/illumina/final_raw_reads/"$F1"_merged_R2.fastq.gz) > /dev/null 2>&1 
else
(mv $path/"$ncbi"data/illumina/raw_reads/"$F1"_*R1*.gz $path/"$ncbi"data/illumina/final_raw_reads/"$F1"_original_R1.fastq.gz) > /dev/null 2>&1 
(mv $path/"$ncbi"data/illumina/raw_reads/"$F1"_*R2*.gz $path/"$ncbi"data/illumina/final_raw_reads/"$F1"_original_R2.fastq.gz) > /dev/null 2>&1 
fi

###############################################################################
if [ ! -d results ]; then
mkdir results
fi

if [ ! -d results/01_raw_read_count ]; then
mkdir results/01_raw_read_count
fi

if [ ! -d results/01_raw_read_count/$F1 ]; then
mkdir results/01_raw_read_count/$F1
fi

#------------------------------------------------------------------------------

echo "running.... raw_read_count for..." $F1 
V1=$(zcat $path/"$ncbi"data/illumina/final_raw_reads/"$F1"_*R1.fastq.gz | echo $((`wc -l`/4))) # total number of F reads
V2=$(zcat $path/"$ncbi"data/illumina/final_raw_reads/"$F1"_*R2.fastq.gz | echo $((`wc -l`/4))) # total number of R reads

zcat $path/"$ncbi"data/illumina/final_raw_reads/"$F1"_*R1.fastq.gz | awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' > results/01_raw_read_count/$F1/$F1.1_average_read_length.tmp.statistics.tab
V3=$(cat results/01_raw_read_count/$F1/$F1.1_average_read_length.tmp.statistics.tab | awk '{$3=$1*$2; C+=$3; B+=$2; ARL= C / B} END {print ARL}')  #ARL
rm results/01_raw_read_count/$F1/$F1.1_average_read_length.tmp.statistics.tab

zcat $path/"$ncbi"data/illumina/final_raw_reads/"$F1"_*R2.fastq.gz | awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' > results/01_raw_read_count/$F1/$F1.1_average_read_length.tmp.statistics.tab
V4=$(cat results/01_raw_read_count/$F1/$F1.1_average_read_length.tmp.statistics.tab | awk '{$3=$1*$2; C+=$3; B+=$2; ARL= C / B} END {print ARL}') #ARL
rm results/01_raw_read_count/$F1/$F1.1_average_read_length.tmp.statistics.tab

V5=$( awk "BEGIN {print ($V3+$V4)/2; exit}") # average of ARL
#V6=$( awk "BEGIN {print (($V1+$V2)*$V5)/5000000; exit}") # Dont remember what is this 
echo $F1 $V1 $V2 $V5 > results/01_raw_read_count/$F1/$F1.1_raw_read_count.statistics.tab

(rm results/01_raw_read_count/raw_read_count.statistics.tab) > /dev/null 2>&1 
echo $F1 $V1 $V2 $V5 >> results/01_raw_read_count/raw_read_count.statistics.tab

if [ ! "$V1" == "$V2" ] ; then
echo "Warnig: $F1 $V1 not_equalt_to $V2"
echo "$F1 $V1 not_equalt_to $V2"  warning.txt
fi

if [ ! -f results/01_raw_read_count/$F1/$F1.1_raw_read_count.statistics.tab ]; then
echo "failed! raw read count for $F1" 
echo "failed! raw read count for $F1" >> results/failed.txt
fi


echo "finished... raw_read_count for..." $F1 
done

###############################################################################

echo "completed.... step-1 raw_read_count ------------------------------------"

###############################################################################

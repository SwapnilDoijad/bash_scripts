#!/bin/bash
###############################################################################
#2 filter_read
###############################################################################
echo "started...... step-2 filter_reads --------------------------------------"
###############################################################################
echo "provide list file (for e.g. all)"
ls list.*.txt | sed 's/ /\n/g'
read l
list=$(echo "list.$l.txt")
echo "raw data path? (for e.g. /media/network/reads_database/p_Entb_Germany)"
ls /media/network/reads_database/
read path_tmp1
path=$(echo "/media/network/reads_database/$path_tmp1")
(mkdir $path/results) > /dev/null 2>&1
(mkdir $path/results/02_filtered_reads) > /dev/null 2>&1
# -----------------------------------------------------------------------------
echo "Do you want to sample the reads,type y else press ENTER to continue"
read answer
if [ "$answer" == "y" ]; then
    echo "what is the expected size of the genome? (for eg. 5000000)"
    read genome_size
    echo "How much genome coverage needed? (for eg. 60)"
    read user_coverage
    for F1 in $(cat $list); do
        echo "running.... step-2 filter_reads for..." $F1 
        (mkdir $path/results/02_filtered_reads/$F1)> /dev/null 2>&1
        (mkdir $path/results/02_filtered_reads/$F1/raw_reads_sampled)> /dev/null 2>&1
        #read sampling ---------------------------------------------------------------------

        ARL=$(awk '{print $4}' results/01_raw_read_count/$F1/"$F1".1_raw_read_count.statistics.tab | awk -F'.' '{print $1}')
        reads_for_desired_coverage=$(( $genome_size * $user_coverage / $ARL ))
        /home/swapnil/pipeline/tools/seqtk/seqtk sample -s100 $path/data/illumina/final_raw_reads/"$F1"_*R1.fastq.gz $reads_for_desired_coverage > $path/results/02_filtered_reads/$F1/raw_reads_sampled/"$F1"_sampled_R1.fastq
        /home/swapnil/pipeline/tools/seqtk/seqtk sample -s100 $path/data/illumina/final_raw_reads/"$F1"_*R2.fastq.gz $reads_for_desired_coverage > $path/results/02_filtered_reads/$F1/raw_reads_sampled/"$F1"_sampled_R2.fastq
        gzip -c $path/results/02_filtered_reads/$F1/raw_reads_sampled/"$F1"_sampled_R1.fastq > $path/results/02_filtered_reads/$F1/raw_reads_sampled/"$F1"_sampled_R1.fastq.gz
        gzip -c $path/results/02_filtered_reads/$F1/raw_reads_sampled/"$F1"_sampled_R2.fastq > $path/results/02_filtered_reads/$F1/raw_reads_sampled/"$F1"_sampled_R2.fastq.gz

        rm $path/results/02_filtered_reads/$F1/raw_reads_sampled/"$F1"_sampled_R1.fastq 
        rm $path/results/02_filtered_reads/$F1/raw_reads_sampled/"$F1"_sampled_R2.fastq 

        ## read sampling: Filter reads by TRIMMOMATIC------------------------------------
        java -jar /home/swapnil/tools/trimmomatic/trimmomatic-0.36.jar PE -threads 8 -phred33 -quiet \
        $path/results/02_filtered_reads/$F1/raw_reads_sampled/"$F1"_sampled_R1.fastq.gz \
        $path/results/02_filtered_reads/$F1/raw_reads_sampled/"$F1"_sampled_R1.fastq.gz \
        $path/results/02_filtered_reads/$F1/"$F1"_R1_001.filtered_paired.fastq.gz \
        $path/results/02_filtered_reads/$F1/"$F1"_R1_001.filtered_unpaired.fastq.gz \
        $path/results/02_filtered_reads/$F1/"$F1"_R2_001.filtered_paired.fastq.gz \
        $path/results/02_filtered_reads/$F1/"$F1"_R2_001.filtered_unpaired.fastq.gz \
        ILLUMINACLIP:/home/swapnil/tools/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 LEADING:6 TRAILING:6 SLIDINGWINDOW:4:20 MINLEN:36
        #### read sampling: FIlter reads by TRIMMOMATIC -------------------------------

        echo "finished... step-2 filter_reads for..." $F1
        done

        #------------------------------------------------------------------------------

        else

        for F1 in $(cat $list); do
        echo "running.... step-2 filter_reads for..." $F1 

        if [ ! -d $path/results/02_filtered_reads/$F1 ]; then
        mkdir $path/results/02_filtered_reads/$F1
        mkdir $path/results/02_filtered_reads/$F1/raw_reads_sampled
        fi

        java -jar /home/swapnil/tools/trimmomatic/trimmomatic-0.36.jar PE -threads 8 -phred33 -quiet \
        $path/data/illumina/final_raw_reads/"$F1"_*R1.fastq.gz \
        $path/data/illumina/final_raw_reads/"$F1"_*R2.fastq.gz \
        $path/results/02_filtered_reads/$F1/"$F1"_R1_001.filtered_paired.fastq.gz \
        $path/results/02_filtered_reads/$F1/"$F1"_R1_001.filtered_unpaired.fastq.gz \
        $path/results/02_filtered_reads/$F1/"$F1"_R2_001.filtered_paired.fastq.gz \
        $path/results/02_filtered_reads/$F1/"$F1"_R2_001.filtered_unpaired.fastq.gz \
        ILLUMINACLIP:/home/swapnil/tools/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 LEADING:6 TRAILING:6 SLIDINGWINDOW:4:20 MINLEN:36
        echo "finished... step-2 filter_reads for..." $F1

    done
fi

#FILE-CHECK--------------------------------------------------------------------

for F1 in $(cat $list); do

if [ ! -f $path/results/02_filtered_reads/$F1/"$F1"_R1_001.filtered_paired.fastq.gz ]; then
echo "filtering reads for $F1 failed" 
echo "filtering reads for $F1 failed" >> $path/results/failed_list.txt
fi

if [ ! -f $path/results/02_filtered_reads/$F1/"$F1"_R2_001.filtered_paired.fastq.gz ]; then
echo "filtering reads for $F1 failed" 
echo "filtering reads for $F1 failed" >> $path/results/failed_list.txt
fi

done

###############################################################################

echo "completed.... step-2 filter_reads --------------------------------------"

###############################################################################

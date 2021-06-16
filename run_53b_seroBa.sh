
#!/bin/bash
###############################################################################
#53b SeroBa
###############################################################################
echo "Hi swapnil, script is not working"
exit
    echo "provide list file (for e.g. all)"
    echo "-------------------------------------------------------------------------"
    ls list.*.txt | sed 's/ /\n/g'
    echo "-------------------------------------------------------------------------"
    read l
    list=$(echo "list.$l.txt")

    echo "raw data path? (for e.g. /media/network/reads_database/p_Entb_Germany)"
    read path

    (mkdir $path/results) > /dev/null 2>&1
    (mkdir $path/results/53b_seroBa) > /dev/null 2>&1
    (mkdir $path/results/53b_seroBa/raw_files) > /dev/null 2>&1
###############################################################################

for F1 in $(cat $list); do
    rename 's/_001.filtered_paired.fastq/\.fq/'  $path/results/02_filtered_reads/$F1/"$F1"_R1_001.filtered_paired.fastq.gz
    rename 's/_001.filtered_paired.fastq/\.fq/'  $path/results/02_filtered_reads/$F1/"$F1"_R2_001.filtered_paired.fastq.gz
    (mkdir $path/results/53b_seroBa/raw_files/$F1) > /dev/null 2>&1
    docker run --rm -it -v /home/ubuntu/data:/data sangerpathogens/seroba seroba runSerotyping seroba/database \
    $path/results/02_filtered_reads/$F1/"$F1"_R1.fq.gz \
    $path/results/02_filtered_reads/$F1/"$F1"_R2.fq.gz \
    $path/results/53b_seroBa/raw_files/$F1
    rename 's/\.fq/_001.filtered_paired.fastq/'  $path/results/02_filtered_reads/$F1/"$F1"_R1.fq.gz
    rename 's/\.fq/_001.filtered_paired.fastq/'  $path/results/02_filtered_reads/$F1/"$F1"_R2.fq.gz
done

###############################################################################
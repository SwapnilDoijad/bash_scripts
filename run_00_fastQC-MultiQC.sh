## Note: fastqc also work with nanopore
###############################################################################
## run fastq and multiQC

    echo "provide list file (for e.g. all)"
    echo "---------------------------------------------------------------------"
    ls list.*.txt | sed 's/ /\n/g'
    echo "---------------------------------------------------------------------"
    read l
    list=$(echo "list.$l.txt")

    echo "Full raw read path? (for e.g. /media/swapnil/network/reads_database/p_Entb_Germany)"
    read path

    echo "do you want to run for raw_reads? type y" 
    read raw_reads 

    echo "do you want to run for filtered_reads? type y" 
    read filtered_reads 

###############################################################################

    (mkdir results) > /dev/null 2>&1
    (mkdir results/00_read_QC) > /dev/null 2>&1
    (mkdir results/00_read_QC/raw_files) > /dev/null 2>&1

    ## for raw reads
    if [ "$raw_reads" == "y" ]; then
        for F1 in $(cat $list); do
        if [ ! -f results/00_read_QC/raw_files/$F1/"$F1"_R1_001_fastqc.html ]; then
        echo "running fastqc for $F1"
        (mkdir results/00_read_QC/raw_files/$F1) > /dev/null 2>&1
            if [ -f $path/data/illumina/raw_reads/"$F1"_R1_001.fastq.gz ] ; then 
                (/home/swapnil/tools/FastQC/fastqc -o results/00_read_QC/raw_files/$F1 -f fastq -q -t 8 $path/data/illumina/raw_reads/"$F1"_*R1*.gz) > /dev/null 2>&1
                (/home/swapnil/tools/FastQC/fastqc -o results/00_read_QC/raw_files/$F1 -f fastq -q -t 8 $path/data/illumina/raw_reads/"$F1"_*R2*.gz) > /dev/null 2>&1
                else
                (/home/swapnil/tools/FastQC/fastqc -o results/00_read_QC/raw_files/$F1 -f fastq -q -t 8 $path/data/illumina/final_reads/"$F1"_*R1*.gz ) > /dev/null 2>&1
                (/home/swapnil/tools/FastQC/fastqc -o results/00_read_QC/raw_files/$F1 -f fastq -q -t 8 $path/data/illumina/final_reads/"$F1"_*R2*.gz ) > /dev/null 2>&1
            fi
            else
            echo "raw-read fastqc report already generated for $F1"
        fi
        done

        source /home/swapnil/miniconda3/etc/profile.d/conda.sh
        conda activate myenv
        echo "running multiqc"
        multiqc -q -o results/00_read_QC/ results/00_read_QC/raw_files/
        conda deactivate
    fi

    ## for filtered reads
    if [ "$filtered_reads" == "y" ]; then
        (mkdir $path/results/02_filtered_reads/raw_files) > /dev/null 2>&1
        for F1 in $(cat $list); do
        (mkdir $path/results/02_filtered_reads/raw_files/$F1) > /dev/null 2>&1
        if [ ! -f $path/results/02_filtered_reads/$F1/"$F1"_R1_001_fastqc.html ]; then
        echo "running fastqc for $F1"
            (/home/swapnil/tools/FastQC/fastqc -o $path/results/02_filtered_reads/raw_files/$F1 -f fastq -q -t 8 $path/results/02_filtered_reads/$F1/"$F1"_R1_001.filtered_paired.fastq.gz) > /dev/null 2>&1
            (/home/swapnil/tools/FastQC/fastqc -o $path/results/02_filtered_reads/raw_files/$F1 -f fastq -q -t 8 $path/results/02_filtered_reads/$F1/"$F1"_R2_001.filtered_paired.fastq.gz) > /dev/null 2>&1
            else
            echo "filtered-read fastqc report already generated for $F1"
        fi
        done

        source /home/swapnil/miniconda3/etc/profile.d/conda.sh
        conda activate myenv
        echo "running multiqc"
        multiqc -q -o $path/results/02_filtered_reads/ $path/results/02_filtered_reads/raw_files/
        conda deactivate
    fi

###############################################################################
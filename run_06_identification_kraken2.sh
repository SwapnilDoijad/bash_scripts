###############################################################################
## preliminary file preparations

    echo "provide list file (for e.g. all)"
    echo "---------------------------------------------------------------------"
    ls list.*.txt | sed 's/ /\n/g'
    echo "---------------------------------------------------------------------"
    read l
    list=$(echo "list.$l.txt")

    echo "Classification based on fasta (f) or raw reads (r) or filted-reads (s)"
    read input_file_type
    sfx=$(echo _"$input_file_type")
    if [  "$input_file_type" == "r" ] || [  "$input_file_type" == "s" ] ; then 
        echo "raw data path? (for e.g. /media/swapnil/network/reads_database/p_Entb_Germany)"
        read path
    fi

    (mkdir results/06_identification_kraken2"$sfx" ) > /dev/null 2>&1 
    (mkdir results/06_identification_kraken2"$sfx"/raw_files ) > /dev/null 2>&1
    path1=results/06_identification_kraken2"$sfx"
    path2=results/06_identification_kraken2"$sfx"/raw_files
###############################################################################
## run kraken2
    if [  "$input_file_type" == "f" ] ; then 
        for F1 in $(cat $list ); do
            echo "running kraken2 (fasta) for $F1"
            ( mkdir results/06_identification_kraken2"$sfx"/raw_files/$F1 ) > /dev/null 2>&1  
            ( kraken2 --db /media/swapnil/share/databases/minikraken2_v1_8GB --threads 8 --report $path2/$F1/report.txt --use-names results/04_assembly/all_fasta/$F1.fasta ) > /dev/null 2>&1 
            #( kraken2 --db /media/swapnil/share/databases/minikraken2_v1_8GB --threads 8 --unclassified-out $path2/$F1/unclassified.txt  --classified-out $path2/$F1/classified.txt --report $path2/$F1/report.txt --use-names --output $path2/$F1/results.txt results/04_assembly/all_fasta/$F1.fasta ) > /dev/null 2>&1 
            #gzip -c $path2/$F1/"unclassified#".fastq > $path2/$F1/"unclassified#".fastq.gz
            #gzip -c $path2/$F1/"classified#".fastq > $path2/$F1/"classified#".fastq.gz
            A1=$(awk -F'\t' ' $4=="S" {print $0}' $path2/$F1/report.txt | head -1 | sed 's/\t                /\t/g')
            echo $F1 $A1 >> $path1/report.csv

        done
        
        elif [ "$input_file_type" == "r" ] ; then
        for F1 in $(cat $list ); do
            echo "running kraken2 (raw_reads) for $F1"
            ( mkdir results/06_identification_kraken2"$sfx"/raw_files/$F1 ) > /dev/null 2>&1 
            ( kraken2 --db /media/swapnil/share/databases/minikraken2_v1_8GB --threads 8 --report $path2/$F1/report.txt --use-names --paired $path/data/illumina/final_reads/"$F1"_*R1*.gz $path/data/illumina/final_reads/"$F1"_*R2*.gz ) > /dev/null 2>&1
            #( kraken2 --db /media/swapnil/share/databases/minikraken2_v1_8GB --threads 8 --unclassified-out $path2/$F1/"unclassified#".fastq --classified-out $path2/$F1/"$F1"_"classified#".fastq --report $path2/$F1/report.txt --use-names --output $path2/$F1/results.txt --paired $path/data/illumina/final_reads/"$F1"_*R1*.gz $path/data/illumina/final_reads/"$F1"_*R2*.gz ) > /dev/null 2>&1 
            #gzip -c $path2/$F1/"unclassified#".fastq > $path2/$F1/"unclassified#".fastq.gz
            #gzip -c $path2/$F1/"classified#".fastq > $path2/$F1/"classified#".fastq.gz
            A1=$(awk -F'\t' ' $4=="S" {print $0}' $path2/$F1/report.txt | head -1 | sed 's/\t                /\t/g')
            echo $F1 $A1 >> $path1/report.csv
        done
        
        elif [  "$input_file_type" == "s" ] ; then
        for F1 in $(cat $list ); do
            echo "running kraken2 (filtered_reads) for $F1"
            ( mkdir results/06_identification_kraken2"$sfx"/raw_files/$F1 ) > /dev/null 2>&1 
            ( kraken2 --db /media/swapnil/share/databases/minikraken2_v1_8GB --threads 8 --report $path2/$F1/report.txt --use-names --paired $path/results/02_filtered_reads/$F1/"$F1"_R1_001.filtered_paired.fastq.gz $path/results/02_filtered_reads/$F1/"$F1"_R2_001.filtered_paired.fastq.gz ) > /dev/null 2>&1 
            #( kraken2 --db /media/swapnil/share/databases/minikraken2_v1_8GB --threads 8 --unclassified-out $path2/$F1/"unclassified#".fastq  --classified-out $path2/$F1/"$F1"_"classified#".fastq --report $path2/$F1/report.txt --use-names --output $path2/$F1/results.txt --paired $path/results/02_filtered_reads/$F1/"$F1"_R1_001.filtered_paired.fastq.gz $path/results/02_filtered_reads/$F1/"$F1"_R2_001.filtered_paired.fastq.gz ) > /dev/null 2>&1 
            #gzip -c $path2/$F1/"unclassified#".fastq > $path2/$F1/"unclassified#".fastq.gz
            #gzip -c $path2/$F1/"classified#".fastq > $path2/$F1/"classified#".fastq.gz
            A1=$(awk -F'\t' ' $4=="S" {print $0}' $path2/$F1/report.txt | head -1 | sed 's/\t                /\t/g')
            echo $F1 $A1 >> $path1/report.csv
        done
    fi
###############################################################################
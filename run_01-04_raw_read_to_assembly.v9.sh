#!/bin/bash
###############################################################################
## raw read count, filter reads, filter read count and assemble
## Warning: Do not save orignial four files of NextSeq  and a single merged file *R1*.fastq.gz in one folder, else the read count will be performed for both
#NanoStat --fastq /media/network/reads_data"$ncbi"base/p_Entb_Germany/data"$ncbi"/nanopore/final_reads/Survcare220_BC12.fastq.gz --outdir /media/swapnil/Swapnil_Sony_HD-E1/data"$ncbi"/Genome_data"$ncbi"base/p_Eb_two_Germany/results/01_raw_read_count/00_nanopore -p Survcare220_BC12 -n Survcare220_BC12
###############################################################################
## preliminary file preparations

    echo "provide list file (for e.g. all)"
    echo "---------------------------------------------------------------------"
    ls list.*.txt | sed 's/ /\n/g' | awk -F'.' '{print $2}'
    echo "---------------------------------------------------------------------"
    read l
    list=$(echo "list.$l.txt")

    echo "raw data"$ncbi" provide full path? (for e.g. /media/swapnil/network/reads_database/p_Entb_Germany/)"
    echo "---------------------------------------------------------------------"
    read path

    if [ -d "$path"data/ncbi ] ; then
    echo "do you want to analyze the NCBI data"$ncbi"? type y for yes" 
    read ncbi_tmp
        if [ "$ncbi_tmp" == "y" ] ; then
            ncbi=$(echo "/ncbi")
            else
            ncbi=$(echo "")
            echo "Do you want to sample the reads,type y else press ENTER to continue"
            read answer
        fi
    fi

    echo "do you want to sample the data? type y" 
    read sampling 
        if [ "$sampling" == "y" ] ; then 
            echo "what is the expected size of the genome? (for eg. 5000000)"
            read genome_size
            echo "How much genome coverage needed? (for eg. 60)"
            read user_coverage
            else
            genome_size=$(echo "5000000")
        fi

    #echo "Do you want to run quast for all fasta? answer yes or PRESS ENTER to skip" 
    #read quast

    (mkdir "$path"data"$ncbi") > /dev/null 2>&1 
    (mkdir "$path"data"$ncbi"/illumina) > /dev/null 2>&1 
    (mkdir "$path"data"$ncbi"/illumina/final_reads) > /dev/null 2>&1

echo "###############################################################################"


## 1 raw read count step
    echo "started...... step-1 raw_read_count -----------------------------------------"
    for F1 in $(cat $list);do
        echo "running $F1"
        ## check if already done
        if [ -f results/01_raw_read_count/raw_files/$F1/$F1.01_raw_read_count.statistics.tab ] ; then
            echo "raw_read_count for $F1 already finished" 
            else
            ## echo "$F1 is not calculated"
            ## short-read count only if short-reads are present
            #echo ""$path"data"$ncbi"/illumina/raw_reads/"$F1"_*R1*.gz"
            total_number_of_reads=$( ls "$path"data"$ncbi"/illumina/raw_reads/"$F1"_*R1*.gz 2>/dev/null | wc -l )
            #echo "$total_number_of_reads total_number_of_reads"
            #if [ -f "$path"data"$ncbi"/illumina/raw_reads/"$F1"_*R1*.gz ] && [ -f "$path"data"$ncbi"/illumina/raw_reads/"$F1"_*R2*.gz ] ; then ##dont know why did not work
            if [ "$total_number_of_reads" -gt 0 ] ; then
                echo "Identify the raw reads of Miseq/NextSeq and place them in final_reads"
                V1=$(ls "$path"data"$ncbi"/illumina/raw_reads/"$F1"*.gz | wc -l )
                #echo "$V1 V1"
                if (( $V1 > 2 )); then
                    cat "$path"data"$ncbi"/illumina/raw_reads/"$F1"_*R1*.gz > "$path"data"$ncbi"/illumina/final_reads/"$F1"_merged_R1.fastq.gz #) > /dev/null 2>&1 
                    cat "$path"data"$ncbi"/illumina/raw_reads/"$F1"_*R2*.gz > "$path"data"$ncbi"/illumina/final_reads/"$F1"_merged_R2.fastq.gz #) > /dev/null 2>&1 
                    #( rm "$path"data"$ncbi"/illumina/raw_reads/"$F1"_*R1*.gz ) > /dev/null 2>&1 
                    #( rm "$path"data"$ncbi"/illumina/raw_reads/"$F1"_*R2*.gz ) > /dev/null 2>&1 
                    else
                    (mv "$path"data"$ncbi"/illumina/raw_reads/"$F1"_*R1*.gz "$path"data"$ncbi"/illumina/final_reads/"$F1"_original_R1.fastq.gz) > /dev/null 2>&1 
                    (mv "$path"data"$ncbi"/illumina/raw_reads/"$F1"_*R2*.gz "$path"data"$ncbi"/illumina/final_reads/"$F1"_original_R2.fastq.gz) > /dev/null 2>&1 
                fi
                else
                echo "$F1 raw reads did not located"
            fi

            if [ -f "$path"data"$ncbi"/illumina/final_reads/"$F1"_*R1*.gz ] && [ -f "$path"data"$ncbi"/illumina/final_reads/"$F1"_*R2*.gz ] ; then
                ## echo "Okay, $F1 is already in final_reads/"
                (mkdir results) > /dev/null 2>&1 
                (mkdir results/01_raw_read_count) > /dev/null 2>&1
                (mkdir results/01_raw_read_count/raw_files) > /dev/null 2>&1 
                (mkdir results/01_raw_read_count/raw_files/$F1) > /dev/null 2>&1 
                #------------------------------------------------------------------------------

                echo "running.... raw_read_count for..." $F1 
                V1=$(zcat "$path"data"$ncbi"/illumina/final_reads/"$F1"_*R1*.fastq.gz | echo $((`wc -l`/4))) # total number of F reads
                V2=$(zcat "$path"data"$ncbi"/illumina/final_reads/"$F1"_*R2*.fastq.gz | echo $((`wc -l`/4))) # total number of R reads

                V3=$(zcat "$path"data"$ncbi"/illumina/final_reads/"$F1"_*R1*.fastq.gz | awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' | awk '{$3=$1*$2; C+=$3; B+=$2; ARL= C / B} END {print ARL}')  #ARL
                V4=$(zcat "$path"data"$ncbi"/illumina/final_reads/"$F1"_*R2*.fastq.gz | awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' | awk '{$3=$1*$2; C+=$3; B+=$2; ARL= C / B} END {print ARL}') #ARL
                total_raw_reads=$(( $V1 + $V2 ))
                #echo total_raw_reads $total_raw_reads

                V5=$( awk "BEGIN {print ($V3+$V4)/2; exit}" | awk -F'.' '{print $1}') # average of ARL 
                raw_read_coverage=$(($total_raw_reads*$V5/$genome_size))
                echo $F1 $V1 $V2 $total_raw_reads $V5 $raw_read_coverage  > results/01_raw_read_count/raw_files/$F1/$F1.01_raw_read_count.statistics.tab

                ## write stat to common stat-result file
                F1_already_present=$(grep "$F1" results/01_raw_read_count/raw_files/$F1/$F1.01_raw_read_count.statistics.tab)
                if [ ! -z "$F1_already_present" ]; then #if stat-data already present
                    echo $F1 $V1 $V2 $total_raw_reads $V5 $raw_read_coverage >> results/01_raw_read_count/raw_read_count.statistics.tab
                fi

                if [ ! "$V1" == "$V2" ] ; then
                    echo "Warnig: $F1 $V1 not_equalt_to $V2"
                    echo "$F1 $V1 not_equalt_to $V2" >> 01_raw_read_warning.txt
                    else
                    :
                fi

                if [ ! -f results/01_raw_read_count/raw_files/$F1/$F1.01_raw_read_count.statistics.tab ]; then
                    echo "failed! raw read count for $F1" 
                    echo "failed! raw read count for $F1" >> results/failed.txt
                    else
                    :
                fi

                ## raw-read count/stat for long-read
                
                elif [ -f "$path"data"$ncbi"/long_read/raw_reads/$F1.fastq.gz ] ; then
                echo "running long read counts for $F1"
                (mkdir results/01_raw_read_count/raw_files/$F1) > /dev/null 2>&1 
                NanoStat --fastq "$path"data"$ncbi"/long_read/raw_reads/"$F1".fastq.gz -o results/01_raw_read_count/raw_files/$F1 -n $F1.01_raw_read_count.long_read.statistics.tab
            fi
        fi

        ## a step for - if raw reads avilable for both short-reads and long-read
        ## short-read, by default will get counted above, while long-read get skipped, therefore this step is necessary
        if [ ! -f results/01_raw_read_count/raw_files/$F1/$F1.01_raw_read_count.statistics.tab ] ; then
            if [ -f "$path"data"$ncbi"/long_read/raw_reads/$F1.fastq.gz ] ; then
                echo "running long read counts for $F1"
                (mkdir results/01_raw_read_count/raw_files/$F1) > /dev/null 2>&1
                NanoStat --fastq "$path"data"$ncbi"/long_read/raw_reads/"$F1".fastq.gz -o results/01_raw_read_count/raw_files/$F1 -n $F1.01_raw_read_count.long_read.statistics.tab
            fi
        fi
    done
###############################################################################
## 2 filter read step
    echo "started...... step-2 filter and count filtered reads ------------------------"
    (mkdir "$path"results) > /dev/null 2>&1
    (mkdir "$path"results/02_filtered_reads) > /dev/null 2>&1
    # -----------------------------------------------------------------------------
    ## as NCBI is already sampled, skip this sampling step for ncbi data & run only for user raw read
    ## the parameters for the raw-reads are in step-1
    if [ "$ncbi_tmp" == "ncbi" ] ; then
        ## only trimming for NCBI data
        for F1 in $(cat $list); do
            if [ -f "$path"data"$ncbi"/illumina/final_reads/"$F1"_*R1*.gz ] && [ -f "$path"data"$ncbi"/illumina/final_reads/"$F1"_*R2*.gz ] ; then
                echo "running.... step-2 filter_reads for..." $F1 
                (mkdir "$path"results/02_filtered_reads/$F1)> /dev/null 2>&1
                (mkdir "$path"results/02_filtered_reads/$F1/raw_reads_sampled)> /dev/null 2>&1
                java -jar /home/swapnil/tools/trimmomatic/trimmomatic-0.36.jar PE -threads 8 -phred33 -quiet \
                "$path"data"$ncbi"/illumina/final_reads/"$F1"_*R1*.gz \
                "$path"data"$ncbi"/illumina/final_reads/"$F1"_*R2*.gz \
                "$path"results/02_filtered_reads/$F1/"$F1"_R1_001.filtered_paired.fastq.gz \
                "$path"results/02_filtered_reads/$F1/"$F1"_R1_001.filtered_unpaired.fastq.gz \
                "$path"results/02_filtered_reads/$F1/"$F1"_R2_001.filtered_paired.fastq.gz \
                "$path"results/02_filtered_reads/$F1/"$F1"_R2_001.filtered_unpaired.fastq.gz \
                ILLUMINACLIP:/home/swapnil/tools/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 LEADING:6 TRAILING:6 SLIDINGWINDOW:4:20 MINLEN:36
                echo "finished... step-2 filter_reads for..." $F1
            fi
        done

        elif [ "$ncbi_tmp" != "ncbi" ] ; then
            for F1 in $(cat $list); do
                if [ -f "$path"results/02_filtered_reads/"$F1"/"$F1"_R1_001.filtered_paired.fastq.gz ] && [ -f "$path"results/02_filtered_reads/"$F1"/"$F1"_R2_001.filtered_paired.fastq.gz ] ; then
                    echo "filtering and filtered read count for $F1 is already finished"
                    else
                    echo "running sampling and raw read count for $F1"
                    if [ -f "$path"data"$ncbi"/illumina/final_reads/"$F1"_*R1*.gz ] && [ -f "$path"data"$ncbi"/illumina/final_reads/"$F1"_*R2*.gz ] ; then
                        echo "running.... step-2 filter_reads for..." $F1 
                        (mkdir "$path"results/02_filtered_reads/$F1)> /dev/null 2>&1
                        (mkdir "$path"results/02_filtered_reads/$F1/raw_reads_sampled)> /dev/null 2>&1
                        ## read sampling ---------------------------------------------------------------------
                        if [ "$sampling" == "y" ] ; then
                            ## check if sampling already done
                            if [ -f "$path"results/02_filtered_reads/"$F1"/raw_reads_sampled/"$F1"_sampled_R1.fastq.gz ] &&  [ -f "$path"results/02_filtered_reads/"$F1"/raw_reads_sampled/"$F1"_sampled_R2.fastq.gz ] ; then 
                                echo "$F1 already sampled"
                                else
                                ARL=$(awk '{print $5}' results/01_raw_read_count/raw_files/$F1/"$F1".01_raw_read_count.statistics.tab | awk '{printf "%.0f\n", $1}' ) 
                                reads_for_desired_coverage=$(( $genome_size * $user_coverage / $ARL ))
                                /home/swapnil/tools/seqtk/seqtk sample -s100 "$path"data"$ncbi"/illumina/final_reads/"$F1"_*R1*.fastq.gz $reads_for_desired_coverage > "$path"results/02_filtered_reads/$F1/raw_reads_sampled/"$F1"_sampled_R1.fastq
                                /home/swapnil/tools/seqtk/seqtk sample -s100 "$path"data"$ncbi"/illumina/final_reads/"$F1"_*R2*.fastq.gz $reads_for_desired_coverage > "$path"results/02_filtered_reads/$F1/raw_reads_sampled/"$F1"_sampled_R2.fastq
                                echo "compressing $F1"
                                gzip -c "$path"results/02_filtered_reads/$F1/raw_reads_sampled/"$F1"_sampled_R1.fastq > "$path"results/02_filtered_reads/$F1/raw_reads_sampled/"$F1"_sampled_R1.fastq.gz
                                gzip -c "$path"results/02_filtered_reads/$F1/raw_reads_sampled/"$F1"_sampled_R2.fastq > "$path"results/02_filtered_reads/$F1/raw_reads_sampled/"$F1"_sampled_R2.fastq.gz

                                rm "$path"results/02_filtered_reads/$F1/raw_reads_sampled/"$F1"_sampled_R1.fastq 
                                rm "$path"results/02_filtered_reads/$F1/raw_reads_sampled/"$F1"_sampled_R2.fastq 
                            fi

                            ## read sampling: Filter reads by TRIMMOMATIC------------------------------------
                            ## check if filtering already done
                            if [ -f "$path"results/02_filtered_reads/"$F1"/"$F1"_R1_001.filtered_paired.fastq.gz ] && [ -f "$path"results/02_filtered_reads/"$F1"/"$F1"_R2_001.filtered_paired.fastq.gz ] ; then
                                echo "reads for $F1 already filtered"
                                else
                                echo "filtering reads for $F1"
                                java -jar /home/swapnil/tools/trimmomatic/trimmomatic-0.36.jar PE -threads 8 -phred33 -quiet \
                                "$path"results/02_filtered_reads/$F1/raw_reads_sampled/"$F1"_*R1.fastq.gz \
                                "$path"results/02_filtered_reads/$F1/raw_reads_sampled/"$F1"_*R1.fastq.gz \
                                "$path"results/02_filtered_reads/$F1/"$F1"_R1_001.filtered_paired.fastq.gz \
                                "$path"results/02_filtered_reads/$F1/"$F1"_R1_001.filtered_unpaired.fastq.gz \
                                "$path"results/02_filtered_reads/$F1/"$F1"_R2_001.filtered_paired.fastq.gz \
                                "$path"results/02_filtered_reads/$F1/"$F1"_R2_001.filtered_unpaired.fastq.gz \
                                ILLUMINACLIP:/home/swapnil/tools/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 LEADING:6 TRAILING:6 SLIDINGWINDOW:4:20 MINLEN:36
                            fi
                            #### read sampling: FIlter reads by TRIMMOMATIC -------------------------------

                            ## if no sampling the go directly filtering raw reads
                            elif [ "$sampling" != "y" ] ; then
                                ## read sampling: Filter reads by TRIMMOMATIC------------------------------------
                                ## check if filtering already done
                                if [ -f "$path"results/02_filtered_reads/"$F1"/"$F1"_R1_001.filtered_paired.fastq.gz ] && [ -f "$path"results/02_filtered_reads/"$F1"/"$F1"_R2_001.filtered_paired.fastq.gz ] ; then
                                    echo "reads for $F1 already filtered"
                                    else
                                    echo "filtering reads for $F1"
                                    java -jar /home/swapnil/tools/trimmomatic/trimmomatic-0.36.jar PE -threads 8 -phred33 -quiet \
                                    "$path"data"$ncbi"/illumina/final_reads/"$F1"_*R1*.fastq.gz \
                                    "$path"data"$ncbi"/illumina/final_reads/"$F1"_*R2*.fastq.gz \
                                    "$path"results/02_filtered_reads/$F1/"$F1"_R1_001.filtered_paired.fastq.gz \
                                    "$path"results/02_filtered_reads/$F1/"$F1"_R1_001.filtered_unpaired.fastq.gz \
                                    "$path"results/02_filtered_reads/$F1/"$F1"_R2_001.filtered_paired.fastq.gz \
                                    "$path"results/02_filtered_reads/$F1/"$F1"_R2_001.filtered_unpaired.fastq.gz \
                                    ILLUMINACLIP:/home/swapnil/tools/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 LEADING:6 TRAILING:6 SLIDINGWINDOW:4:20 MINLEN:36
                                fi
                            #### read sampling: FIlter reads by TRIMMOMATIC -------------------------------
                        fi
                        echo "finished... step-2 filter_reads for..." $F1
                    fi
                fi
            done
    fi
###############################################################################
## 3 filtered read count step
    echo "started...... step-3 filtered read count-------------------------------------"
    (mkdir results/03_filtered_reads_count) > /dev/null 2>&1 
    ##-----------------------------------------------------------------------------
    for F1 in $(cat $list); do
        ## check if already finished
        if [ -f results/03_filtered_reads_count/"$F1"/tmp/"$F1".3_filtered_read_count.statistics.tab ] ; then
            echo "filtered read count already finished for... $F1"
            else
            ## short-read count only if short-reads are present
            if [ -f "$path"data"$ncbi"/illumina/final_reads/"$F1"_*R1*.gz ] && [ -f "$path"data"$ncbi"/illumina/final_reads/"$F1"_*R2*.gz ] ; then
                
                if [ ! -f "$path"results/03_filtered_reads_count/$F1/"$F1".3_filtered_read_count.statistics.tab ]; then

                    echo "running.... step-3 filtered_read_count for..." $F1
                    (mkdir results/03_filtered_reads_count/$F1) > /dev/null 2>&1 
                    (mkdir results/03_filtered_reads_count/$F1/tmp) > /dev/null 2>&1 

                    V1=$(zcat "$path"results/02_filtered_reads/$F1/"$F1"_R1_001.filtered_paired.fastq.gz | echo $((`wc -l`/4)))
                    V2=$(zcat "$path"results/02_filtered_reads/$F1/"$F1"_R2_001.filtered_paired.fastq.gz | echo $((`wc -l`/4)))

                    zcat "$path"results/02_filtered_reads/$F1/"$F1"_R1_001.filtered_paired.fastq.gz | awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' > results/03_filtered_reads_count/$F1/tmp/$F1.3_filtered_read_count.tmp.statistics.tab
                    V3=$(cat results/03_filtered_reads_count/$F1/tmp/$F1.3_filtered_read_count.tmp.statistics.tab | awk '{$3=$1*$2; C+=$3; B+=$2; ARL= C / B} END {print ARL}') 
                    rm results/03_filtered_reads_count/$F1/tmp/$F1.3_filtered_read_count.tmp.statistics.tab

                    zcat "$path"results/02_filtered_reads/$F1/"$F1"_R2_001.filtered_paired.fastq.gz | awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' > results/03_filtered_reads_count/$F1/tmp/$F1.3_filtered_read_count.tmp.statistics.tab
                    V4=$(cat results/03_filtered_reads_count/$F1/tmp/$F1.3_filtered_read_count.tmp.statistics.tab | awk '{$3=$1*$2; C+=$3; B+=$2; ARL= C / B} END {print ARL}') 
                    rm results/03_filtered_reads_count/$F1/tmp/$F1.3_filtered_read_count.tmp.statistics.tab

                    V5=$(awk "BEGIN {print ($V3+$V4)/2; exit}" | awk -F'.' '{print $1}' ) #ARL
                    total_filtered_reads=$(( $V1 + $V2 ))
                    filtered_read_coverage=$(( $total_filtered_reads*$V5/$genome_size ))

                    echo $F1 $V1 $V2 $V5 $filtered_read_coverage > results/03_filtered_reads_count/$F1/tmp/$F1.3_filtered_read_count.statistics.tab

                    (V6=$(grep "$F1" results/03_filtered_reads_count/filtered_read_count.statistics.tab)) > /dev/null 2>&1 
                    if [ "$V6" !=  "$F1" ]; then
                        echo $F1 $V1 $V2 $V5 $filtered_read_coverage >> results/03_filtered_reads_count/filtered_read_count.statistics.tab
                    fi

                    if [ ! -f results/03_filtered_reads_count/$F1/tmp/$F1.3_filtered_read_count.statistics.tab ]; then
                        echo "failed! $F1 filtered read count " 
                        echo "failed! $F1 filtered read count " >> results/failed.txt
                    fi

                    echo "finished... step-3 filtered_read_count for..." $F1
                    else
                    echo "Reads already filtered for $F1 "
                fi
            fi
        fi
    done
###############################################################################
## 4 assembly step
    echo "started...... step-4 assembly -----------------------------------------------"
    (mkdir results/04_assembly) > /dev/null 2>&1 
    (mkdir results/04_assembly/raw_files) > /dev/null 2>&1 
    (mkdir results/04_assembly/all_fasta) > /dev/null 2>&1 
    
    # check short and long read, if both available
    for F1 in $(cat $list); do
        if [ -f "$path"results/02_filtered_reads/$F1/"$F1"_R1_001.filtered_paired.fastq.gz ] && [ -f "$path"data"$ncbi"/long_read/raw_reads/$F1.fastq.gz ] ; then
            ##  unicycler hybrid assembly
                if [ ! -f results/04_assembly/unicycler/$F1/$F1.fasta ]; then
                    echo "running hybrid assembly for $F1"
                    WorDir=$(echo $PWD)
                    /home/swapnil/tools/Unicycler_v0.4.8/unicycler-runner.py \
                    -1 "$path"results/02_filtered_reads/$F1/"$F1"_R1_001.filtered_paired.fastq.gz \
                    -2 "$path"results/02_filtered_reads/$F1/"$F1"_R2_001.filtered_paired.fastq.gz \
                    -l "$path"data"$ncbi"/long_read/raw_reads/"$F1".fastq.gz \
                    -o $WorDir/results/04_assembly/unicycler/$F1 \
                    --spades_path /home/swapnil/tools/SPAdes-3.12.0/bin/spades.py \
                    --pilon_path /home/swapnil/tools/pilon-1.22.jar \
                    --min_fasta_length 500 \
                    --keep 0 \
                    --verbosity 0
                    cp results/04_assembly/unicycler/$F1/assembly.fasta results/04_assembly/unicycler/$F1/$F1.joined.fasta
                    sed -i 's/>.*/NNNNNNNNNN/g' results/04_assembly/unicycler/$F1/$F1.joined.fasta
                    sed -i "1i "'>'$F1"" results/04_assembly/unicycler/$F1/$F1.joined.fasta
                    sed -i "\$aNNNNNNNNNN" results/04_assembly/unicycler/$F1/$F1.joined.fasta
                    cp results/04_assembly/unicycler/$F1/$F1.joined.fasta results/04_assembly/all_fasta/$F1.fasta
                    else
                    echo "hybrid assembly already finished for $F1"
                fi
            elif [ -f "$path"data"$ncbi"/long_read/raw_reads/$F1.fastq.gz ] ; then
            ##  unicycler long-read assembly
                if [ ! -f results/04_assembly/unicycler/$F1/$F1.fasta ]; then
                    echo "running long-read assembly for $F1"
                    WorDir=$(echo $PWD)
                    /home/swapnil/tools/Unicycler_v0.4.8/unicycler-runner.py \
                    -l "$path"data"$ncbi"/long_read/raw_reads/"$F1".fastq.gz \
                    -o $WorDir/results/04_assembly/unicycler/$F1 \
                    --spades_path /home/swapnil/tools/SPAdes-3.12.0/bin/spades.py \
                    --pilon_path /home/swapnil/tools/pilon-1.22.jar \
                    --min_fasta_length 500 \
                    --keep 0 \
                    --verbosity 0
                    cp results/04_assembly/unicycler/$F1/assembly.fasta results/04_assembly/unicycler/$F1/$F1.joined.fasta
                    sed -i 's/>.*/NNNNNNNNNN/g' results/04_assembly/unicycler/$F1/$F1.joined.fasta
                    sed -i "1i "'>'$F1"" results/04_assembly/unicycler/$F1/$F1.joined.fasta
                    sed -i "\$aNNNNNNNNNN" results/04_assembly/unicycler/$F1/$F1.joined.fasta
                    cp results/04_assembly/unicycler/$F1/$F1.joined.fasta results/04_assembly/all_fasta/$F1.fasta
                else
                echo "long-read assembly already finished for $F1"
                fi
            elif [ -f "$path"results/02_filtered_reads/$F1/"$F1"_R1_001.filtered_paired.fastq.gz ] && [ -f "$path"results/02_filtered_reads/$F1/"$F1"_R2_001.filtered_paired.fastq.gz ] ; then
            ##  SPADES assembly by unicycler
                if [ ! -f results/04_assembly/raw_files/$F1/$F1.fasta ]; then
                    echo "running SPADES assembly by unicycler for $F1"
                    WorDir=$(echo $PWD)
                    /home/swapnil/tools/Unicycler_v0.4.8/unicycler-runner.py \
                    -1 "$path"results/02_filtered_reads/$F1/"$F1"_R1_001.filtered_paired.fastq.gz \
                    -2 "$path"results/02_filtered_reads/$F1/"$F1"_R2_001.filtered_paired.fastq.gz \
                    -o $WorDir/results/04_assembly/raw_files/$F1 \
                    --spades_path /home/swapnil/tools/SPAdes-3.15.1-Linux/bin/spades.py \
                    --pilon_path /home/swapnil/tools/pilon-1.22.jar \
                    --min_fasta_length 500 \
                    --keep 0 
                    # --verbosity 0 ## just put off this option for moment
                    (mkdir $WorDir/results/04_assembly/raw_files/$F1/tmp) > /dev/null 2>&1
                    rm -rf results/04_assembly/raw_files/$F1/configs
                    rm -rf results/04_assembly/raw_files/$F1/corrected
                    rm -rf results/04_assembly/raw_files/$F1/pilon_polish
                    cp results/04_assembly/raw_files/$F1/assembly.fasta results/04_assembly/raw_files/$F1/contigs.fasta ## this is to get similar output file of contigs.fasta
                    cp results/04_assembly/raw_files/$F1/contigs.fasta results/04_assembly/raw_files/$F1/$F1.fasta
                    total_contigs=$(grep -c ">" results/04_assembly/raw_files/$F1/$F1.fasta) ## total number of contigs in
                    genome_length=$(grep ">" results/04_assembly/raw_files/$F1/$F1.fasta | awk -F'=' '{print $2}' | awk '{print $1}' | awk '{sum+=$1;} END{print sum;}')
                    echo $total_contigs $genome_length > results/04_assembly/raw_files/$F1/tmp/$F1.statistics.tab
                    cp results/04_assembly/raw_files/$F1/assembly.fasta results/04_assembly/raw_files/$F1/$F1.joined.fasta
                    sed -i 's/>.*/NNNNNNNNNN/g' results/04_assembly/raw_files/$F1/$F1.joined.fasta
                    sed -i "1i "'>'$F1"" results/04_assembly/raw_files/$F1/$F1.joined.fasta
                    sed -i "\$aNNNNNNNNNN" results/04_assembly/raw_files/$F1/$F1.joined.fasta
                    cp results/04_assembly/raw_files/$F1/$F1.joined.fasta results/04_assembly/all_fasta/$F1.fasta
                    ###############################################################################
                    ## Quality control by assembly-stat and quast
                    (assembly-stats -t  results/04_assembly/raw_files/$F1/assembly.fasta >  results/04_assembly/raw_files/$F1/$F1.raw-assembly.stat.tab) > /dev/null 2>&1
                    (/home/swapnil/tools/quast-5.0.2/quast.py  results/04_assembly/raw_files/$F1/$F1.fasta --output-dir  results/04_assembly/raw_files/$F1/quast_raw_fasta) > /dev/null 2>&1
                    bioawk -c fastx '{ print $name, length($seq) }' <results/04_assembly/raw_files/$F1/contigs.fasta > results/04_assembly/raw_files/$F1/contigs_size.csv
                    #echo "finished... step-4 assembly for..." $F1
                    else
                    echo "Assembly already finished for $F1"
                fi
            else
            echo "$F1 Long-read or Long-short-read absent PLEASE CHECK ---------------"
            echo "$F1 Long-read or Long-short-read absent PLEASE CHECK ---------------" > 04_assembly/unicycler.failed.list
        fi
        ## reformat fasta 
        fasta_formatter -i results/04_assembly/all_fasta/$F1.fasta -w 60 -o results/04_assembly/all_fasta/$F1.fasta.tmp
        mv results/04_assembly/all_fasta/$F1.fasta.tmp results/04_assembly/all_fasta/$F1.fasta
        ## this script also can be used to reformat fasta, check input 
        # (python /home/swapnil/tools/read-cleaning-format-conversion/KSU_bioinfo_lab/fasta-o-matic/fasta_o_matic.py -o results/04_assembly/all_fasta/ -f results/04_assembly/all_fasta/$F1.fasta)> /dev/null 2>&1
    done
    ###############################################################################
    ## All-together quast
    if [ "$quast" == "yes" ]; then
        if [ ! -f results/00_ref/ref.*.fasta ] && [ ! -f results/00_ref/ref.*.gff ]; then
            echo "ref.*.fasta or ref.*.gff is ABSENT, please add them in results/00_ref/ and then press enter"
            read -p "Press enter to continue"
        fi

        (mkdir results/04_assembly/00_icarus) > /dev/null 2>&1
        (mkdir results/04_assembly/00_icarus/contigs) > /dev/null 2>&1

        for F1 in $(cat $list); do
            cp  results/04_assembly/raw_files/$F1/$F1.fasta results/04_assembly/00_icarus/contigs/
        done
        (quast.py results/04_assembly/00_icarus/contigs/*.fasta -R results/00_ref/ref.*.gbk -G results/00_ref/ref.*.gff --output-dir results/04_assembly/00_icarus/) > /dev/null 2>&1
    fi
###############################################################################
## create stat file
    (rm results/04_assembly/raw_files/contigs.stat.csv) > /dev/null 2>&1
    (rm results/04_assembly/contigs.stat.xlsx) > /dev/null 2>&1
    (touch results/04_assembly/raw_files/contigs.stat.csv) > /dev/null 2>&1
    echo "running stat"
    for F1 in $(cat $list); do
    V1=$(cat results/04_assembly/raw_files/$F1/tmp/$F1.statistics.tab)
    (V2=$(cat results/04_assembly/raw_files/contigs.stat.csv | grep "$F1" )) > /dev/null 2>&1
            if [ -z "$V2" ] ; then
                    echo $F1 $V1 
                    echo $F1 $V1 >> results/04_assembly/raw_files/contigs.stat.csv
            else
            echo "$F1 already in contigs.stat.csv"
            fi
    done
    sed -i '1 i\isolate >500bp_contigs assembly_size' results/04_assembly/raw_files/contigs.stat.csv
    sed -i 's/\s\+/,/g' results/04_assembly/raw_files/contigs.stat.csv
    ssconvert results/04_assembly/raw_files/contigs.stat.csv results/04_assembly/contigs.stat.xlsx
###############################################################################
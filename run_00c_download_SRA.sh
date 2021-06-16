#!/bin/bash
###############################################################################
#00 download and process SRA
###############################################################################
echo "started...... step-00 download and process SRA -------------------------"
###############################################################################
## inital input and directory creation step
    echo “Update: What is the name of the genus species?”
    read name
    echo "want to store the data on external drive, then provide path? (for e.g. /media/swapnil/network/reads_database/p_Enterobacter_NCBI)"
	read path_tmp
	if [ ! -z $path_tmp ] ; then
		path=$(echo "$path_tmp/")	
        else
        path=$(echo "$path_tmp")
	fi
    #------------------------------------------------------------------------------
    ## genomic platform details
    LS=$(echo "GENOMIC") #LibrarySource
    LL=$(echo "SINGLE") #LibraryLayout
    Plf=$(echo "PACBIO_SMRT") #Platform
    Plfm=$(echo "$Plf" | tr '[:upper:]' '[:lower:]')
    #------------------------------------------------------------------------------
    (mkdir data) > /dev/null 2>&1
    (mkdir data/ncbi) > /dev/null 2>&1
    (mkdir data/ncbi/$Plfm ) > /dev/null 2>&1
    (mkdir data/ncbi/$Plfm/raw_reads ) > /dev/null 2>&1

##############################################################################
## Download SRA
    for SRABiosample in $(cat data/ncbi/list."$name".SRA-BioSample.to-update.$LS.$LL.$Plf.txt); do
    file_count=$(find "$path"data/ncbi/illumina -name "$SRABiosample"_*.gz 2>/dev/null | wc -l)
    if [ $file_count -gt 0 ] ; then
    echo "$SRABiosample already finished"
    else
        SRABiosample_run=$(awk -F'\t' '{print $2}' data/ncbi/metadata/tmp/"$SRABiosample"_SRA.metadata.$LS.$LL.$Plf.tab)
        SRABiosample_size=$(awk -F'\t' '{print $5}' data/ncbi/metadata/tmp/"$SRABiosample"_SRA.metadata.$LS.$LL.$Plf.tab)
        echo "downloading $SRABiosample ($SRABiosample_size Mb) from SRA"
        (/home/swapnil/tools/sratoolkit.2.8.1-centos_linux64/bin/prefetch $SRABiosample_run) > /dev/null 2>&1 ;
        echo "converting SRA to fastq $SRABiosample"
        (/home/swapnil/tools/sratoolkit.2.8.1-centos_linux64/bin/fastq-dump  --outdir /home/swapnil/ncbi/$Plfm/raw_reads/ --split-files /home/swapnil/ncbi/public/sra/$SRABiosample_run.sra) > /dev/null 2>&1 ;
        if [ -f /home/swapnil/ncbi/$Plfm/raw_reads/"$SRABiosample_run"_1.fastq ] && [ -f /home/swapnil/ncbi/$Plfm/raw_reads/"$SRABiosample_run"_2.fastq ] ; then
        (mkdir data/ncbi/$Plfm/raw_reads/$SRABiosample) > /dev/null 2>&1
        cp /home/swapnil/ncbi/$Plfm/raw_reads/"$SRABiosample_run"_1.fastq data/ncbi/$Plfm/raw_reads/$SRABiosample/"$SRABiosample"_1.fastq
        rm /home/swapnil/ncbi/$Plfm/raw_reads/"$SRABiosample_run"_1.fastq
        cp /home/swapnil/ncbi/$Plfm/raw_reads/"$SRABiosample_run"_2.fastq data/ncbi/$Plfm/raw_reads/$SRABiosample/"$SRABiosample"_2.fastq
        rm /home/swapnil/ncbi/$Plfm/raw_reads/"$SRABiosample_run"_2.fastq
        rm /home/swapnil/ncbi/public/sra/$SRABiosample_run.sra
    ##############################################################################
    ## preliminary count for the average read length and sample
        echo "sampling and processing $SRABiosample"
        (mkdir data/ncbi/$Plfm/raw_reads/$SRABiosample/tmp) > /dev/null 2>&1
        C1=$(cat data/ncbi/$Plfm/raw_reads/$SRABiosample/"$SRABiosample"_1.fastq | echo $((`wc -l`/4)))
        cat data/ncbi/$Plfm/raw_reads/$SRABiosample/"$SRABiosample"_1.fastq | awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' > data/ncbi/$Plfm/raw_reads/$SRABiosample/tmp/$SRABiosample.read-statistics.1.tmp
        cat data/ncbi/$Plfm/raw_reads/$SRABiosample/tmp/$SRABiosample.read-statistics.1.tmp | awk '{$3=$1*$2; C+=$3; B+=$2; ARL= C / B} END {print ARL}' > data/ncbi/$Plfm/raw_reads/$SRABiosample/tmp/$SRABiosample.read-statistics.2.tmp
        V1=$(awk -F'.' '{print $1}' data/ncbi/$Plfm/raw_reads/$SRABiosample/tmp/$SRABiosample.read-statistics.2.tmp)
        V2=$(( 125000000 / $V1 )) # calculate the number of reads required for 50x coverage for 5Mb genome (25x coverage reads F and 25X coverage reads R) (this will be applied to sample the F or R reads sufficient for 50x)
        /home/swapnil/tools/seqtk/seqtk sample -s100 data/ncbi/$Plfm/raw_reads/$SRABiosample/"$SRABiosample"_1.fastq $V2 > data/ncbi/$Plfm/raw_reads/$SRABiosample/"$SRABiosample"_sampled_R1.fastq
        /home/swapnil/tools/seqtk/seqtk sample -s100 data/ncbi/$Plfm/raw_reads/$SRABiosample/"$SRABiosample"_2.fastq $V2 > data/ncbi/$Plfm/raw_reads/$SRABiosample/"$SRABiosample"_sampled_R2.fastq
        gzip -f data/ncbi/$Plfm/raw_reads/$SRABiosample/"$SRABiosample"_sampled_R1.fastq > data/ncbi/$Plfm/raw_reads/$SRABiosample/"$SRABiosample"_sampled_R1.fastq.gz
        gzip -f data/ncbi/$Plfm/raw_reads/$SRABiosample/"$SRABiosample"_sampled_R2.fastq > data/ncbi/$Plfm/raw_reads/$SRABiosample/"$SRABiosample"_sampled_R2.fastq.gz
        echo $V1 $C1 > data/ncbi/$Plfm/raw_reads/$SRABiosample/tmp/$SRABiosample.read-statistics.tmp
    ##############################################################################
    ## for SAMPLED data: count for the average read length
        C2=$(zcat data/ncbi/$Plfm/raw_reads/$SRABiosample/"$SRABiosample"_sampled_R1.fastq.gz | echo $((`wc -l`/4)))
        zcat data/ncbi/$Plfm/raw_reads/$SRABiosample/"$SRABiosample"_sampled_R1.fastq.gz | awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' > data/ncbi/$Plfm/raw_reads/$SRABiosample/tmp/$SRABiosample.sampled.read-statistics.1.tmp
        cat data/ncbi/$Plfm/raw_reads/$SRABiosample/tmp/$SRABiosample.read-statistics.1.tmp | awk '{$3=$1*$2; C+=$3; B+=$2; ARL= C / B} END {print ARL}' > data/ncbi/$Plfm/raw_reads/$SRABiosample/tmp/$SRABiosample.sampled.read-statistics.2.tmp
        V3=$(awk -F'.' '{print $1}' data/ncbi/$Plfm/raw_reads/$SRABiosample/tmp/$SRABiosample.sampled.read-statistics.2.tmp)
        echo $V1 $C2 > data/ncbi/$Plfm/raw_reads/$SRABiosample/tmp/$SRABiosample.sampled.read-statistics.tmp

        rm data/ncbi/$Plfm/raw_reads/$SRABiosample/"$SRABiosample"_1.fastq
        rm data/ncbi/$Plfm/raw_reads/$SRABiosample/"$SRABiosample"_2.fastq

        A1=$(awk '{print $1}' data/ncbi/$Plfm/raw_reads/$SRABiosample/tmp/$SRABiosample.read-statistics.tmp) #read-length
        A2=$(awk '{print $2}' data/ncbi/$Plfm/raw_reads/$SRABiosample/tmp/$SRABiosample.read-statistics.tmp) #total-reads
        total_coverage=$(( 2 * $A1 * $A2 / 5000000 )) #total-read-coverage ##need to work on this
        A3=$(awk '{print $1}' data/ncbi/$Plfm/raw_reads/$SRABiosample/tmp/$SRABiosample.sampled.read-statistics.tmp) #sampled-read-length
        A4=$(awk '{print $2}' data/ncbi/$Plfm/raw_reads/$SRABiosample/tmp/$SRABiosample.sampled.read-statistics.tmp) #sampled-total-reads
        sampled_coverage=$(( 2 * $A3 * $A4 / 5000000 )) ##need to work on 
        echo $SRABiosample $A1 $A2 $total_coverage $A3 $A4 $sampled_coverage >> data/ncbi/$Plfm/raw_reads/stat.tmp
        sort -u data/ncbi/$Plfm/raw_reads/stat.tmp > data/ncbi/$Plfm/raw_reads/stat.txt
        #------------------------------------------------------------------------------
        ## move the files to network or external drive
            if [ ! -z $path ] ; then
            echo moving files to $path
            (mkdir "$path"data ) > /dev/null 2>&1
            (mkdir "$path"data/ncbi ) > /dev/null 2>&1
            (mkdir "$path"data/ncbi/$Plfm ) > /dev/null 2>&1
            (mkdir "$path"data/ncbi/$Plfm/raw_reads ) > /dev/null 2>&1
            cp -rf data/ncbi/$Plfm/raw_reads/$SRABiosample/*.gz "$path"data/ncbi/$Plfm/raw_reads/
            rm -rf data/ncbi/$Plfm/raw_reads/$SRABiosample
            fi

        else
        echo "failed $SRABiosample" >> list.failed.SRA-pair-absent.txt
        fi
    fi
    done
###############################################################################
echo "Completed...... step-00 download and process SRA -----------------------"
###############################################################################
## Download SRA long reads
    (mkdir data) > /dev/null 2>&1
    (mkdir data/ncbi) > /dev/null 2>&1
    (mkdir data/ncbi/long_read) > /dev/null 2>&1
    (mkdir data/ncbi/long_read/raw_reads) > /dev/null 2>&1
    for SRABiosample in $(cat data/ncbi/list."$name".SRA-BioSample.to-update.$LS.SINGLE.LONG_READ.txt); do
        file_count=$(find "$path"data/ncbi/long_read/ -name "$SRABiosample"*.gz 2>/dev/null | wc -l)
        if [ $file_count -gt 0 ] ; then
        echo "$SRABiosample already finished"
        else
            SRABiosample_run=$(awk -F'\t' '{print $2}' data/ncbi/metadata/tmp/"$SRABiosample"_SRA.metadata.$LS.SINGLE.LONG_READ.tab)
            SRABiosample_size=$(awk -F'\t' '{print $5}' data/ncbi/metadata/tmp/"$SRABiosample"_SRA.metadata.$LS.SINGLE.LONG_READ.tab)
            echo "downloading $SRABiosample ($SRABiosample_size Mb) from SRA"
            (/home/swapnil/tools/sratoolkit.2.8.1-centos_linux64/bin/prefetch $SRABiosample_run) > /dev/null 2>&1 ;
            echo "converting SRA ($SRABiosample_run) to fastq (biosample $SRABiosample)"
            (/home/swapnil/tools/sratoolkit.2.8.1-centos_linux64/bin/fastq-dump --outdir /home/swapnil/ncbi/long_read/raw_reads/ /home/swapnil/ncbi/public/sra/$SRABiosample_run.sra) > /dev/null 2>&1 ;
            if [ -f /home/swapnil/ncbi/long_read/raw_reads/"$SRABiosample_run".fastq ] ; then
            (mkdir data/ncbi/long_read/raw_reads/$SRABiosample) > /dev/null 2>&1
            cp /home/swapnil/ncbi/long_read/raw_reads/"$SRABiosample_run".fastq data/ncbi/long_read/raw_reads/$SRABiosample/"$SRABiosample".fastq
            gzip data/ncbi/long_read/raw_reads/$SRABiosample/"$SRABiosample".fastq
            rm /home/swapnil/ncbi/long_read/raw_reads/"$SRABiosample_run".fastq
            rm /home/swapnil/ncbi/public/sra/$SRABiosample_run.sra
        ##############################################################################
            ## move the files to netwrok or external drive
                if [ ! -z $path ] ; then
                    echo moving files to $path
                    (mkdir "$path"data ) > /dev/null 2>&1
                    (mkdir "$path"data/ncbi ) > /dev/null 2>&1
                    (mkdir "$path"data/ncbi/long_read ) > /dev/null 2>&1
                    (mkdir "$path"data/ncbi/long_read/raw_reads ) > /dev/null 2>&1
                    cp -rf data/ncbi/long_read/raw_reads/$SRABiosample/*.gz "$path"data/ncbi/long_read/raw_reads/
                    rm -rf data/ncbi/long_read/raw_reads/$SRABiosample
                fi
            else
            echo "failed $SRABiosample" >> list.failed.SRA-pair-absent.txt
            fi
        fi
    done
###############################################################################
echo "Completed...... step-00 download and process SRA -----------------------"
###############################################################################
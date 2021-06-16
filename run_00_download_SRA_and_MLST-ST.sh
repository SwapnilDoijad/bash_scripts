#!/bin/bash
###############################################################################
#00 download and process SRA

###############################################################################

echo "started...... step-00 download SRA and MLST ST -------------------------"

###############################################################################

if [ ! -d /media/network/a_baumannii/data ]; then
mkdir /media/network/a_baumannii/data
fi

if [ ! -d /media/network/a_baumannii/data/ncbi-fastq ]; then
mkdir /media/network/a_baumannii/data/ncbi-fastq
fi

if [ ! -d /media/network/a_baumannii/data/ncbi-fastq/raw_reads ]; then
mkdir /media/network/a_baumannii/data/ncbi-fastq/raw_reads
fi

if [ ! -d /media/network/a_baumannii/data/ncbi-fastq/mlst ]; then
mkdir /media/network/a_baumannii/data/ncbi-fastq/mlst
fi

for F1 in $(cat list.SRA.txt);do

## Download

echo "downloding $F1"
(/home/swapnil/pipeline/tools/sratoolkit/bin/prefetch $F1) > /dev/null 2>&1 ;

## Conversion

echo "converting SRA to fastq $F1"
(/home/swapnil/pipeline/tools/sratoolkit/bin/fastq-dump --outdir /home/swapnil/ncbi/fastq/ --split-files /home/swapnil/ncbi/public/sra/$F1.sra) > /dev/null 2>&1 ;

## Check for fastq pair

if [ -f /home/swapnil/ncbi/fastq/"$F1"_1.fastq ] && [ -f /home/swapnil/ncbi/fastq/"$F1"_2.fastq ] ; then

cp /home/swapnil/ncbi/fastq/"$F1"_1.fastq /media/network/a_baumannii/data/ncbi-fastq/raw_reads ;
(gzip /media/network/a_baumannii/data/ncbi-fastq/raw_reads/"$F1"_1.fastq) > /dev/null 2>&1 ;
rm /home/swapnil/ncbi/fastq/"$F1"_1.fastq ;
cp /home/swapnil/ncbi/fastq/"$F1"_2.fastq /media/network/a_baumannii/data/ncbi-fastq/raw_reads/ ;
(gzip /media/network/a_baumannii/data/ncbi-fastq/raw_reads/"$F1"_2.fastq) > /dev/null 2>&1 ;
rm /home/swapnil/ncbi/fastq/"$F1"_2.fastq ;
rm /home/swapnil/ncbi/public/sra/$F1.sra ;

## Run MLST

echo "running MLST for $F1"

tmp_path=$(echo $PWD) ;

(python /home/swapnil/tools/srst2-master/build/lib.linux-x86_64-2.7/srst2/srst2.py --output /media/network/a_baumannii/data/ncbi-fastq/mlst/$F1/$F1 --input_pe /media/network/a_baumannii/data/ncbi-fastq/raw_reads/"$F1"_1.fastq.gz /media/network/a_baumannii/data/ncbi-fastq/raw_reads/"$F1"_2.fastq.gz --mlst_db /home/swapnil/tools/srst2-master/mlst_profile/A_baumannii/Acinetobacter_baumannii#1.fasta --mlst_definitions /home/swapnil/tools/srst2-master/mlst_profile/A_baumannii/abaumannii.txt --mlst_delimiter '_' --mlst_max_mismatch 3) > /dev/null 2>&1

awk 'FNR==2 {print $0}' /media/network/a_baumannii/data/ncbi-fastq/mlst/$F1/"$F1"_*_results.txt >> /media/network/a_baumannii/data/ncbi-fastq/mlst/0_mlst_all.csv ;

rm -rf /media/network/a_baumannii/data/ncbi-fastq/mlst/$F1

#desired_ST=$(cat data/mlst_st.txt) ;
#echo "desired_ST $desired_ST" ;
#ST_of_strain=$(awk 'FNR==2 {print $2}' data/NCBI-fastq/mlst/$F1/"$F1"_*_results.txt) ;
#echo "ST_of_strain $ST_of_strain"

#	if [ "$desired_ST" -eq "$ST_of_strain" ] ; then
#	mkdir data/NCBI-fastq/mlst_ST_"$desired_ST"/$F1 ;
#	cp data/NCBI-fastq/"$F1"_1.fastq data/NCBI-fastq/mlst_ST_"$desired_ST"/$F1 ;
#	cp data/NCBI-fastq/"$F1"_1.fastq data/NCBI-fastq/mlst_ST_"$desired_ST"/$F1 ;
#	cp data/NCBI-fastq/mlst/$F1/$F1 data/NCBI-fastq/mlst_ST_"$desired_ST"/ ;
#	else
#	fi

#cp data/NCBI-fastq/"$F1"*.fastq /media/network/a_baumannii/data/ncbi-fastq


## ----

echo "Download, conversion and MLST finished for $F1 -------------------------"


else
echo "failed $F1" >> list.failed.SRA-pair-absent.txt
fi

done

###############################################################################

echo "Completed...... step-00 download SRA and MLST ST -----------------------"

###############################################################################


#!/bin/bash
###############################################################################
#4 assembly
###############################################################################
echo "started.... step-4 assembly --------------------------------------------"
###############################################################################
## file and directory preparation
	echo "provide list file (for e.g. all)"
	ls list.*.txt | sed 's/ /\n/g'
	read l
	list=$(echo "list.$l.txt")

	echo "raw data path? (for e.g. /media/network/reads_database/p_Entb_Germany)"
	ls /media/network/reads_database/
	read path_tmp1
	path=$(echo "/media/network/reads_database/$path_tmp1")

	(mkdir results/04_assembly) > /dev/null 2>&1 
	(mkdir results/04_assembly/raw_files) > /dev/null 2>&1 
	(mkdir results/04_assembly/all_fasta) > /dev/null 2>&1 
###############################################################################
for F1 in $(cat $list); do

echo "running.... step-4 assembly for..." $F1

if [ ! -d results/04_assembly/raw_files/$F1 ]; then
mkdir results/04_assembly/raw_files/$F1
fi

/home/swapnil/tools/SPAdes-3.12.0/bin/spades.py \
-1 $path/results/02_filtered_reads/$F1/"$F1"_R1_001.filtered_paired.fastq.gz \
-2 $path/results/02_filtered_reads/$F1/"$F1"_R2_001.filtered_paired.fastq.gz \
--cov-cutoff 5 --phred-offset 33 -o results/04_assembly/raw_files/$F1  > /dev/null 2>&1
(rm -r results/04_assembly/raw_files/$F1/K21) > /dev/null 2>&1
(rm -r results/04_assembly/raw_files/$F1/K33) > /dev/null 2>&1
(rm -r results/04_assembly/raw_files/$F1/K55) > /dev/null 2>&1
(rm -r results/04_assembly/raw_files/$F1/K77) > /dev/null 2>&1
(rm -r results/04_assembly/raw_files/$F1/K99) > /dev/null 2>&1
(rm -r results/04_assembly/raw_files/$F1/K127) > /dev/null 2>&1
(rm -r results/04_assembly/raw_files/$F1/corrected) > /dev/null 2>&1
(rm -r results/04_assembly/raw_files/$F1/configs) > /dev/null 2>&1
(rm -r results/04_assembly/raw_files/$F1/misc) > /dev/null 2>&1
(rm -r results/04_assembly/raw_files/$F1/tmp) > /dev/null 2>&1

if [ ! -d results/04_assembly/raw_files/$F1/tmp ]; then
mkdir results/04_assembly/raw_files/$F1/tmp
fi

if [ ! -f results/04_assembly/raw_files/$F1/contigs.fasta ]; then
echo "failed! assembly for $F1" 
echo "failed! assembly for $F1" >> results/failed_list.txt
fi
###############################################################################
cp  results/04_assembly/raw_files/$F1/contigs.fasta  results/04_assembly/raw_files/$F1/$F1.fasta
V1=$(grep -c ">"  results/04_assembly/raw_files/$F1/$F1.fasta)
grep -F ">"  results/04_assembly/raw_files/$F1/$F1.fasta | sed -e 's/_/ /g' | sort -nrk 6 | awk '$6>=5.0 && $4>=500 {print $0}' | sed -s 's/ /_/g' | sed -e 's/>//g' >  results/04_assembly/raw_files/$F1/tmp/$F1.10-500-filtered-contigs.csv
perl /home/swapnil/pipeline/tools/fastagrep.pl -f  results/04_assembly/raw_files/$F1/tmp/$F1.10-500-filtered-contigs.csv  results/04_assembly/raw_files/$F1/$F1.fasta >  results/04_assembly/raw_files/$F1/$F1.contigs-filtered.fasta
V2=$(grep -c ">"  results/04_assembly/raw_files/$F1/$F1.contigs-filtered.fasta)
V3=$(awk -F '_' '{ sum += $6; n++ } END { if (n > 0) print sum / n; }'  results/04_assembly/raw_files/$F1/tmp/$F1.10-500-filtered-contigs.csv)
echo $V1 $V2 $V3 >  results/04_assembly/raw_files/$F1/tmp/$F1.5_filtering_contigs.statistics.tab
cp  results/04_assembly/raw_files/$F1/$F1.contigs-filtered.fasta results/04_assembly/all_fasta/$F1.fasta
sed -i 's/>.*/NNNNNNNNNN/g' results/04_assembly/all_fasta/$F1.fasta
sed -i "\$aNNNNNNNNNN" results/04_assembly/all_fasta/$F1.fasta
sed -i "1i "'>'$F1"" results/04_assembly/all_fasta/$F1.fasta
if [ ! -f  results/04_assembly/raw_files/$F1/tmp/$F1.5_filtering_contigs.statistics.tab ]; then
echo "failed! filtering contigs for $F1" 
echo "failed! filtering contigs for $F1" >> results/failed_list.txt
fi
###############################################################################
## Quality control by assembly-stat and quast
	(assembly-stats -t  results/04_assembly/raw_files/$F1/contigs.fasta >  results/04_assembly/raw_files/$F1/$F1.raw-assembly.stat.tab) > /dev/null 2>&1
	(results/04_assembly/raw_files/$F1/$F1.filtered-assembly.stat.tab) > /dev/null 2>&1
	(assembly-stats -u  results/04_assembly/raw_files/$F1/$F1.contigs-filtered.fasta >>  results/04_assembly/raw_files/$F1/$F1.filtered-assembly.stat.tab) > /dev/null 2>&1
	(quast.py  results/04_assembly/raw_files/$F1/$F1.fasta --output-dir  results/04_assembly/raw_files/$F1/quast_raw_fasta) > /dev/null 2>&1
	(quast.py  results/04_assembly/raw_files/$F1/$F1.contigs-filtered.fasta --output-dir  results/04_assembly/raw_files/$F1/quast_filtered_fasta) > /dev/null 2>&1
	echo "finished... step-4 assembly for..." $F1
	done
###############################################################################
## All-together quast
	echo "Do you want to run quast for all fasta? answer yes or PRESS ENTER to skip" 
	read F2
	if [ "$F2" == "yes" ]; then
		if [ ! -f results/00_ref/ref.*.fasta ] && [ ! -f results/00_ref/ref.*.gff ]; then
		echo "ref.*.fasta or ref.*.gff is ABSENT, please add them in results/00_ref/ and then press enter"
		read -p "Press enter to continue"
		fi

	(mkdir results/04_assembly/00_icarus) > /dev/null 2>&1
	(mkdir results/04_assembly/00_icarus/contigs_filtered_fasta) > /dev/null 2>&1

	for F1 in $(cat $list); do
	cp  results/04_assembly/raw_files/$F1/$F1.contigs-filtered.fasta results/04_assembly/00_icarus/contigs_filtered_fasta/
	done
	(quast.py results/04_assembly/00_icarus/contigs_filtered_fasta/*.fasta -R results/00_ref/ref.*.gbk -G results/00_ref/ref.*.gff --output-dir results/04_assembly/00_icarus/) > /dev/null 2>&1
	fi
###############################################################################
echo "completed.. step-4 assembly -------------------------------------------" 
###############################################################################




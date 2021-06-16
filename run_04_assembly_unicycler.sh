#!/bin/bash
###############################################################################
## long read assembly by unicycler
###############################################################################
echo "started.... step-4c long read assembly UNICYCLER -----------------------"
###############################################################################
(mkdir results/04_assembly) > /dev/null 2>&1
(mkdir results/04_assembly/unicycler) > /dev/null 2>&1
(mkdir results/04_assembly/unicycler/all_fasta) > /dev/null 2>&1

echo "provide list file (for e.g. all)"
read l
list=$(echo "list.$l.txt")

echo "do you want to analyze the NCBI data"$ncbi"? type ncbi for yes" 
read ncbi_tmp
	if [ $ncbi_tmp = "ncbi" ] ; then
		ncbi=$(echo "/$ncbi_tmp")
		else
		ncbi=$(echo "")
	fi

echo "raw data path? (for e.g. /media/network/reads_database/p_Entb_Germany)"
read path
###############################################################################
## unicycler assembly
# check short and long read, if both available
    for F1 in $(cat $list); do
		if [ -f $path/results/02_filtered_reads/$F1/"$F1"_R1_001.filtered_paired.fastq.gz ] && [ -f $path/data"$ncbi"/long_read/raw_reads/"$SRABiosample_run".fastq.gz ] ; then
			echo "running hybrid assembly for $F1"
			WorDir=$(echo $PWD)
			cd /home/swapnil/tools/Unicycler-master
			unicycler \
			-1 $path/results/02_filtered_reads/$F1/"$F1"_R1_001.filtered_paired.fastq.gz \
			-2 $path/results/02_filtered_reads/$F1/"$F1"_R2_001.filtered_paired.fastq.gz \
			-l $path/data"$ncbi"/long_reads/final_reads/"$F1".fastq.gz \
			-o $WorDir/results/04_assembly/unicycler/$F1 \
			--spades_path /home/swapnil/tools/SPAdes-3.12.0/bin/spades.py \
			--pilon_path /home/swapnil/tools/pilon-1.22.jar \
			--keep 0
			cd $WorDir
			cp results/04_assembly/unicycler/$F1/assembly.fasta results/04_assembly/unicycler/$F1/$F1.joined.fasta
			sed -i 's/>.*/NNNNNNNNNN/g' results/04_assembly/unicycler/$F1/$F1.joined.fasta
			sed -i "1i "'>'$F1"" results/04_assembly/unicycler/$F1/$F1.joined.fasta
			sed -i "\$aNNNNNNNNNN" results/04_assembly/unicycler/$F1/$F1.joined.fasta
			cp results/04_assembly/unicycler/$F1/$F1.joined.fasta results/04_assembly/unicycler/all_fasta/$F1.fasta
		elif
			if [ -f $path/data"$ncbi"/long_read/raw_reads/"$SRABiosample_run".fastq.gz ] ; then
			echo "running long-read assembly for $F1"
			WorDir=$(echo $PWD)
			cd /home/swapnil/tools/Unicycler-master
			unicycler \
			-l $path/data"$ncbi"/long_reads/final_reads/"$F1".fastq.gz \
			-o $WorDir/results/04_assembly/unicycler/$F1 \
			--spades_path /home/swapnil/tools/SPAdes-3.12.0/bin/spades.py \
			--pilon_path /home/swapnil/tools/pilon-1.22.jar \
			--keep 0
			cd $WorDir
			cp results/04_assembly/unicycler/$F1/assembly.fasta results/04_assembly/unicycler/$F1/$F1.joined.fasta
			sed -i 's/>.*/NNNNNNNNNN/g' results/04_assembly/unicycler/$F1/$F1.joined.fasta
			sed -i "1i "'>'$F1"" results/04_assembly/unicycler/$F1/$F1.joined.fasta
			sed -i "\$aNNNNNNNNNN" results/04_assembly/unicycler/$F1/$F1.joined.fasta
			cp results/04_assembly/unicycler/$F1/$F1.joined.fasta results/04_assembly/unicycler/all_fasta/$F1.fasta
		else
		echo "Long-read or Long-short-read absent PLEASE CHECK ---------------"
		echo "Long-read or Long-short-read absent PLEASE CHECK ---------------" > 04_assembly/unicycler.failed.list
		done
	fi
###############################################################################
## raw-read stat for PacBio & Nanopore
for F1 in $(cat $list); do
	NanoStat --fastq $path/data"$ncbi"/long_read/final_reads/"$F1".fastq.gz -o results/04_assembly/unicycler/$F1 -n $F1.raw-read.report.txt
done
###############################################################################
## mapping against the reference
echo "Do you want to map the contigs against the reference? yes or PRESS enter to skip the mapping"
read answer
if [ "$answer" == "yes" ] ; then
	if [ -f results/00_ref/ref.*.fasta ]; then
	cp results/00_ref/ref.*.fasta /home/swapnil/tools/mauve/
	else
	echo "missing results/00_ref/ref.*.fasta"
	fi
	for F1 in $(cat $list); do
		echo "mapping $F1 against reference"

		(rm /home/swapnil/tools/mauve/assembly.fasta) > /dev/null 2>&1 
		cp results/04_assembly/unicycler/$F1/assembly.fasta /home/swapnil/tools/mauve/

		WorDir=$(echo $PWD)

		cd /home/swapnil/tools/mauve
		(java -Xmx5000m -cp Mauve.jar org.gel.mauve.contigs.ContigOrderer -output $F1 -ref ref.*.fasta -draft assembly.fasta) > /dev/null 2>&1 

		cd $WorDir

		Var1=$(ls /home/swapnil/tools/mauve/$F1/ | sort -r | head -1)
		cp /home/swapnil/tools/mauve/$F1/$Var1/assembly.fasta  results/04_assembly/unicycler/$F1/"$F1".aligned.fasta
		cp results/04_assembly/unicycler/$F1/"$F1".aligned.fasta results/04_assembly/unicycler/$F1/"$F1".aligned.joined.fasta
		sed -i 's/>.*/NNNNNNNNNN/g' results/04_assembly/unicycler/$F1/"$F1".aligned.joined.fasta
		sed -i "1i "'>'$F1"" results/04_assembly/unicycler/$F1/"$F1".aligned.joined.fasta

		(mkdir results/04_assembly/unicycler/all_fasta) > /dev/null 2>&1

		cp results/04_assembly/unicycler/$F1/"$F1".aligned.joined.fasta results/04_assembly/unicycler/all_fasta/$F1.fasta

		rm -rf /home/swapnil/tools/mauve/ref.$F1.fasta
		rm /home/swapnil/tools/mauve/assembly.fasta

		# stat
		grep ">" results/04_assembly/unicycler/$F1/assembly.fasta | sed 's/>//g' | sed 's/length=//g' | sed 's/depth=//g' > results/04_assembly/unicycler/$F1/stat.csv
		sed -i "1icontig length cov" results/04_assembly/unicycler/$F1/stat.csv
	done
fi
###############################################################################
echo "completed.... step-4 long read UNICYCLER ------------------------------"
###############################################################################
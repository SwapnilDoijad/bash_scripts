#!/bin/bash
###############################################################################
# https://github.com/kblin/ncbi-genome-download

###############################################################################

echo "Downloading genome by ncbi-genome-download------------------------------"

###############################################################################
if [ ! -d data ]; then
mkdir data
fi

if [ ! -d results ]; then
mkdir results
fi

if [ ! -d results/04_assembly ]; then
mkdir results/04_assembly
fi

if [ ! -d results/04_assembly/all_fasta ]; then
mkdir results/04_assembly/all_fasta
fi

#------------------------------------------------------------------------------
## input

echo "assembly-level (all,complete,chromosome,scaffold,contig)"
read assembly_level

echo "genus species ("Serratia marcescens")"
read input_genus

genus=$(echo $input_genus | sed 's/ /_/g')

#------------------------------------------------------------------------------
	## dry rubn to check how many assemblies in TOTAL available to download 
	#ncbi-genome-download -n -g "$input_genus" bacteria

ncbi-genome-download --assembly-level $assembly_level --format fasta -o data/ncbi-genome-download_"$genus"_"$assembly_level" --genus "$input_genus" bacteria

#------------------------------------------------------------------------------

ls data/ncbi-genome-download_"$genus"_"$assembly_level"/refseq/bacteria/ | sed 's/ //g' > list.ncbi-genome-download_"$genus"_"$assembly_level".txt

for F1 in $(cat list.ncbi-genome-download_"$genus"_"$assembly_level".txt); do
gunzip -k data/ncbi-genome-download_"$genus"_"$assembly_level"/refseq/bacteria/"$F1"/*.gz
grep ">" data/ncbi-genome-download_"$genus"_"$assembly_level"/refseq/bacteria/"$F1"/"$F1"*.fna | sed 's/>//g' > data/ncbi-genome-download_"$genus"_"$assembly_level"/refseq/bacteria/"$F1"/list.tmp
awk -F':' '{print $1}' data/ncbi-genome-download_"$genus"_"$assembly_level"/refseq/bacteria/"$F1"/list.tmp | awk -F'.' '{print $1}' > data/ncbi-genome-download_"$genus"_"$assembly_level"/refseq/bacteria/"$F1"/list.2.tmp

perl /home/swapnil/pipeline/tools/split_multifasta.pl --input_file data/ncbi-genome-download_"$genus"_"$assembly_level"/refseq/bacteria/"$F1"/"$F1"_*.fna --output_dir data/ncbi-genome-download_"$genus"_"$assembly_level"/refseq/bacteria/"$F1"/

WrkDir=$(echo $PWD)
cd data/ncbi-genome-download_"$genus"_"$assembly_level"/refseq/bacteria/$F1
	for F3 in $(cat list.2.tmp); do
	V1=$(grep -v ">" "$F3".*.fsa | wc | awk '{print $3-$1}')
	if [ "$V1" -gt 1000000 ] ; then
	mv "$F3".*.fsa "$F3".chr.fsa
	else
	mv "$F3".*.fsa "$F3".plasmid.fsa
	fi
	done
cd $WrkDir

done

###############################################################################

echo "Complete downloading genome by ncbi-genome-download --------------------"

###############################################################################
## arrnage in results/04_assembly/

#for F1 in $(cat list.ncbi-genome-download_"$genus"_"$assembly_level".txt); do
#F2=$(echo $F1 | awk -F'.' '{print $1}' )
#mkdir results/04_assembly/$F2
#mv data/ncbi-genome-download_"$genus"_"$assembly_level"/refseq/bacteria/"$F1"/"$F1"_*.fna results/04_assembly/$F2/$F2.fna
#cp results/04_assembly/$F2/$F2.fna results/04_assembly/all_fasta/$F2.fasta
#done

###############################################################################

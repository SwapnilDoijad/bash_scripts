#!/bin/bash
###############################################################################
#14 WGS-alignment-mugsy


SECONDS=0
###############################################################################

echo "started.... step-14 WGS-alignment-mugsy ---------------------------------"

###############################################################################

if [ ! -d results/14_mugsy_alignment ]; then
mkdir results/14_mugsy_alignment
fi

if [ ! -d results/14_mugsy_alignment/fasta ]; then
mkdir results/14_mugsy_alignment/fasta
fi


for F1 in $(cat list.txt)
do
cp results/04_assembly/all_fasta/$F1.fasta results/14_mugsy_alignment/fasta/
#cp results/0_ref/ref.*.fasta results/14_mugsy_alignment/fasta/
done

#==============================================================================
WorDir=$(echo $PWD)

source /home/swapnil/tools/mugsy/mugsyenv.sh

/home/swapnil/tools/mugsy/mugsy -p alignment --directory "$WorDir"/results/14_mugsy_alignment results/14_mugsy_alignment/fasta/*.fasta > /dev/null 2>&1

rm alignment.mugsy.log

#==============================================================================
## convert mugsy maf format to fasta and extracting snp-sites

V1=$(ls results/14_mugsy_alignment/fasta/*.fasta | wc -l)

awk -f /home/swapnil/pipeline/tools/maf2phy.awk -v n=$V1 results/14_mugsy_alignment/alignment.maf > results/14_mugsy_alignment/alignment.phy

python /home/swapnil/pipeline/tools/convert_sequence_format.mugsy.py

snp-sites -m -o results/14_mugsy_alignment/alignment.snp-sites.fasta results/14_mugsy_alignment/alignment.fasta

#==============================================================================

if [ ! -f results/14_mugsy_alignment/alignment.maf ]; then
echo "failed! mugsy"
echo "failed! mugsy" >> results/failed.txt
fi

###############################################################################

echo "completed...  step-14 WGS-alignment-mugsy -------------------------------"

###############################################################################
duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
###############################################################################

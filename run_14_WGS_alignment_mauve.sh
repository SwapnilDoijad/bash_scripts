#!/bin/bash
###############################################################################
#14 WGS-alignment-mauve


SECONDS=0
###############################################################################

echo "started.... step-14 WGS-alignment-mauve ---------------------------------"

###############################################################################

if [ ! -d results/14_mauve_alignment ]; then
mkdir results/14_mauve_alignment
fi

if [ ! -d results/14_mauve_alignment/fasta ]; then
mkdir results/14_mauve_alignment/fasta
fi


for F1 in $(cat list.txt)
do
cp results/04_assembly/all_fasta/$F1.fasta results/14_mauve_alignment/fasta/
done

#==============================================================================

progressiveMauve --output=results/14_mauve_alignment/out.xmfa results/14_mauve_alignment/fasta/*.fasta

#==============================================================================
## convert mauve xmaf format to fasta and extracting snp-sites

perl /home/swapnil/pipeline/tools/xmfa2fasta.pl --file results/14_mauve_alignment/out.xmfa > results/14_mauve_alignment/out.xmf.fasta

snp-sites -m -o results/14_mauve_alignment/alignment.snp-sites.fasta results/14_mauve_alignment/alignment.fasta

###############################################################################

echo "completed...  step-14 WGS-alignment-mauve -------------------------------"

###############################################################################

## mauve alignmetn can be used for ClonalFrame and gubbins

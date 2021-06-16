#!/bin/bash
###############################################################################
#0 copying_and_converting_data

###############################################################################

echo "started...... step-00 copying_and_converting_data -----------------------"

###############################################################################

echo "reference files should be prefixed with ref (for e.g. ref.file.gb)"

if [ ! -d results ]; then
mkdir results
fi

#------------------------------------------------------------------------------
#processing MAIN reference

if [ ! -d results/00_ref ]; then
mkdir results/00_ref
fi

if [ -s data/references/*.gb ] ; then
rename 's/\.gb/.gbk/' data/references/*.gb
fi

cd data/references
F1=$(ls *.gbk | sed 's/\.gbk//g' | sed 's/ref\.//g')
cd ..
cd ..
#--------------------------------------
# make .gbk file in proper .gbk format

(rm results/00_ref/ref.$F1.gbk)> /dev/null 2>&1

ugene --task=/home/swapnil/pipeline/tools/gbk_to_gbk.uwl --in=data/references/ref."$F1".gbk --out=results/00_ref/ref.$F1.gbk

#--------------------------------------

#python2.7 /home/swapnil/pipeline/tools/faa_extracter_from_gbk.py data/references/ref.$F1.gbk 			#gbk to faa
#mv data/references/$F1.faa results/00_ref/
#sed -i 's/ .*//g' results/00_ref/$F1.faa

perl /home/swapnil/pipeline/tools/GB2Fasta.pl data/references/ref.$F1.gbk results/00_ref/ref.$F1.fasta
echo $F1 | sed -i "1i "'>'$F1"" results/00_ref/ref.$F1.fasta

#------------------------------------------------------------------------------
#for all other ref genome

if [ -s data/references/all_other_ref_genomes/ ] 
then
	(rename 's/\.gb/.gbk/' data/references/all_other_ref_genomes/*.gb) > /dev/null 2>&1
	
	cd data/references/all_other_ref_genomes
	ls *.gbk > all_other_ref_genomes.list.txt
	sed -i 's/\.gbk//g' all_other_ref_genomes.list.txt
	cd ..
	cd ..

	for F2 in $(cat data/references/all_other_ref_genomes/all_other_ref_genomes.list.txt);do
		echo "Started conversion of  $F2"

		ugene --task=/home/swapnil/pipeline/tools/gbk_to_gbk.uwl --in=data/references/all_other_ref_genomes/ref."$F2.gbk" --out=results/00_ref/ref.$F2.gbk

		perl /home/swapnil/pipeline/tools/GB2Fasta.pl data/references/all_other_ref_genomes/ref.$F2.gbk results/00_ref/ref.$F2.fasta
		echo $F2 | sed -i "1i "'>'$F2"" results/00_ref/ref.$F2.fasta
		
		echo "finished conversion of $F2"

	done
fi

#------------------------------------------------------------------------------

#FILE-CHECK--------------------------------------------------------------------

if [ ! -f results/00_ref/ref.$F1.gbk ]; then
echo "failed! copying ref.$F1.gbk" 
echo "failed! copying ref.$F1.gbk" >> results/failed.txt
fi

if [ ! -f results/00_ref/ref.$F1.fasta ]; then
echo "failed! copying ref.$F1.fasta"
echo "failed! copying ref.$F1.fasta" >> results/failed.txt
fi

###############################################################################

echo "completed.... step-00 copying_and_converting_data -----------------------"
###############################################################################

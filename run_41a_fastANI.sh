#!/bin/bash
###############################################################################
#41 fastANI

###############################################################################

echo "started.... step-41 pyANI ----------------------------------------------"

###############################################################################
if [ ! -d results/41a_fastANI ] ; then
mkdir results/41a_fastANI
fi

if [ ! -d results/41a_fastANI/tmp ] ; then
mkdir results/41a_fastANI/tmp
fi

#------------------------------------------------------------------------------

cp list.txt results/41a_fastANI/tmp/list.fastANI.txt

sed -i '/^\s*$/d' results/41a_fastANI/tmp/list.fastANI.txt

sed -i "s/^/results\/04_assembly\/all_fasta\//" results/41a_fastANI/tmp/list.fastANI.txt

sed -i "s/$/\.fasta/" results/41a_fastANI/tmp/list.fastANI.txt

## @SWAPNIL need to work on -matrix option. This option will replace matrix-creation step below
fastANI --ql results/41a_fastANI/tmp/list.fastANI.txt --rl results/41a_fastANI/tmp/list.fastANI.txt -o results/41a_fastANI/tmp/output.csv --matrix

sed -i "s/results\/04_assembly\/all_fasta\///g" results/41a_fastANI/tmp/output.csv

sed -i "s/\.fasta//g" results/41a_fastANI/tmp/output.csv

#------------------------------------------------------------------------------
## matrix-creation step

for F1 in $(cat list.txt);do

	for F2 in $(cat list.txt); do
	awk '$1 == "'$F1'" && $2 == "'$F2'" {print 100-$3}' results/41a_fastANI/tmp/output.csv >> results/41a_fastANI/tmp/$F1.ANI.tmp
	done

sed -i '1 i\'$F1'' results/41a_fastANI/tmp/$F1.ANI.tmp

done

paste results/41a_fastANI/tmp/*.ANI.tmp > results/41a_fastANI/tmp/all.ANI.tmp

awk 'NR==1{print}' results/41a_fastANI/tmp/all.ANI.tmp | /home/swapnil/pipeline/tools/transpose.sh > results/41a_fastANI/tmp/isolate-name.txt
sed -i '1s/^/\n/' results/41a_fastANI/tmp/isolate-name.txt

paste results/41a_fastANI/tmp/isolate-name.txt results/41a_fastANI/tmp/all.ANI.tmp > results/41a_fastANI/all.ANI.tab

Rscript /home/swapnil/pipeline/tools/plot_and_tree.fastANI-distance.r
mv Rplots.pdf results/41a_fastANI/
mv tree.nwk results/41a_fastANI/

#------------------------------------------------------------------------------
# get a representative of a clutster at 95% (-h 5)

echo "running ggrasp.R"

/home/swapnil/tools/GGRaSP-master/ggrasp.R -i results/41a_fastANI/all.ANI.tab -d 100 -h 5 -o results/41a_fastANI/representative_of_a_clusters_at_95_ANI --plottree 

rm Rplots.pdf
###############################################################################

echo "completed.... step-41 pyANI --------------------------------------------"

###############################################################################

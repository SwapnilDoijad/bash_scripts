#!/bin/bash
###############################################################################
#9 WGS-alignment-mugsy

###############################################################################

echo "started.... step-16 bcgTree --------------------------------------------"

###############################################################################
if [ -f results/tmp/faa_list.txt ] ; then
rm results/tmp/faa_list.txt
fi

V1=$(find results/4_annotation/all_faa/ -type f | wc -l); 
i=1
while [ $i -le $V1 ]; 
do 
echo "--proteome bac"$i=; i=$((i+1)); 
done > results/tmp/faa_list.txt


if [ -f results/tmp/all_faa_list.txt ] ; then
rm results/tmp/all_faa_list.txt
fi


for F1 in $(cat list.txt)
do
echo results/4_annotation/all_faa/$F1.faa >> results/tmp/all_faa_list.txt
done

paste -d ""  results/tmp/faa_list.txt results/tmp/all_faa_list.txt > results/tmp/config.txt


if [ ! -d results/16_bcgTree ]; then
mkdir results/16_bcgTree
fi

perl /home/swapnil/pipeline/tools/bcgTree/bin/bcgTree.pl @results/tmp/config.txt --outdir results/16_bcgTree > /dev/null 2>&1


if [ ! -f results/16_bcgTree/RAxML_bipartitions.final ]; then
echo "failed! bcgTree"
echo "failed! bcgTree" >> results/failed_list.txt
fi

###############################################################################

echo "completed.... step-16 bcgTree -------------------------------------------"

###############################################################################

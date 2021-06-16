#!/bin/bash
###############################################################################
# WGS-alignment-mugsy
# add -z option to keep the intermediate files

###############################################################################

echo "started.... badtrip ----------------------------------------------------"

###############################################################################
if [ ! -d results/35_badtrip ] ; then
mkdir results/35_badtrip
fi

if [ ! -d results/35_badtrip/final_data_to_run_badtrip ] ; then
mkdir results/35_badtrip/final_data_to_run_badtrip
fi

if [ ! -d results/35_badtrip/SNV_alignment ] ; then
mkdir results/35_badtrip/SNV_alignment
fi

cp * results/35_badtrip/SNV_alignment

(python3 /home/swapnil/tools/badtrip/format_conversion_scripts/FastaToCounts.py results/35_badtrip/SNV_alignment/*.fasta results/35_badtrip/SNV_alignment/out.txt)> /dev/null 2>&1 ;

cp results/35_badtrip/SNV_alignment/out.txt results/35_badtrip/SNV_alignment/out.tmp ;
sed '1d' results/35_badtrip/SNV_alignment/out.tmp | awk '{$1=$2=""; print $0}' | sed 's/^ *//' | sed 's/,/-/g' | sed 's/ /\t/g' > results/35_badtrip/SNV_alignment/out.2.tmp ;
grep ">" results/35_badtrip/SNV_alignment/*.fasta | sed 's/>//g' > results/35_badtrip/SNV_alignment/list.names.tmp ;

index=1
for V1 in $(cat results/35_badtrip/SNV_alignment/list.names.tmp); do
sed -e "s/$V1/S$index/" -i results/35_badtrip/SNV_alignment/out.2.tmp
        index=$(($index + 1))
done

read -p "place the epidemiological-data.txt and sample-data.txt files in final_data_to_run_badtrip folder and press ENTER " -n1 -s ;

cp results/35_badtrip/SNV_alignment/out.2.tmp results/35_badtrip/final_data_to_run_badtrip/sequence-data.txt

sudo python /home/swapnil/tools/badtrip/scripts/create_BADTRIP_xml.py -a results/35_badtrip/final_data_to_run_badtrip/sequence-data.txt -e results/35_badtrip/final_data_to_run_badtrip/epidemiological-data.txt -s results/35_badtrip/final_data_to_run_badtrip/sample-data.txt -o results/35_badtrip/final_data_to_run_badtrip/for_beast -E True

read -p "If you want to change the .xml file, and press ENTER " -n1 -s ;

sudo java -cp /home/swapnil/tools/badtrip/dist/BADTRIP.v0.1.0/BADTRIP.v0.1.0.src.jar:/home/swapnil/tools/beastV2/lib/beast.jar beast.app.beastapp.BeastMain results/35_badtrip/final_data_to_run_badtrip/for_beast.xml

sudo chown swapnil *.txt
sudo chown swapnil *.log
sudo chown swapnil *.xml

python /home/swapnil/tools/badtrip/scripts/Make_transmission_tree_alternative.py -i results/35_badtrip/final_data_to_run_badtrip/for_beast.trees -o results/35_badtrip/final_data_to_run_badtrip/for_beast.trees.images


###############################################################################

echo "completed.... badtrip --------------------------------------------------"

###############################################################################

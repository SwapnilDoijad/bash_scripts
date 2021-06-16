#!/bin/bash
###############################################################################

###############################################################################

echo "started.... scotti -----------------------------------------------------"

###############################################################################
if [ ! -d results/36_scotti ] ; then
mkdir results/36_scotti
fi

if [ ! -d results/36_scotti/final_data_to_run_scotti ] ; then
mkdir results/36_scotti/final_data_to_run_scotti
fi


read -p "place the sequence.fa dates.csv hosts.csv and hostTimes.csv files in final_data_to_run_scotti folder and press ENTER " -n1 -s ;

sudo python /home/swapnil/tools/scotti/scripts/SCOTTI_generate_xml.py --fasta results/36_scotti/final_data_to_run_scotti/sequence.fa --dates results/36_scotti/final_data_to_run_scotti/dates.csv --hosts results/36_scotti/final_data_to_run_scotti/hosts.csv --hostTimes results/36_scotti/final_data_to_run_scotti/hostTimes.csv --output results/36_scotti/final_data_to_run_scotti/out --maxHosts 40 --numIter 10000 --tracelog 200 --treelog 200 --screenlog 200

sudo java -cp /home/swapnil/tools/scotti/dist/SCOTTI.v1.1.1/SCOTTI.v1.1.1.src.jar:/home/swapnil/tools/beastV2.4.7/lib/beast.jar beast.app.beastapp.BeastMain results/36_scotti/final_data_to_run_scotti/out.xml

python /home/swapnil/tools/scotti/scripts/Make_transmission_tree_alternative.py -i out.trees  -o images


###############################################################################

echo "completed.... scotti ---------------------------------------------------"

###############################################################################


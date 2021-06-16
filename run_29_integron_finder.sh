#!/bin/bash
###############################################################################
#29 integron finder

###############################################################################

echo "started.... integron finder --------------------------------------------"

###############################################################################

if [ ! -d results-YIN_new/29_integron_finder ]; then
mkdir results-YIN_new/29_integron_finder
fi

for F1 in $(cat list-YIN_new1.txt);do

echo "running integron finder for $F1"

(integron_finder results-YIN_new/2_assembly/all_fasta/$F1.fasta --local_max --func_annot --outdir results-YIN_new/29_integron_finder/$F1) > /dev/null 2>&1 

echo "finished integron finder for $F1"
done

###############################################################################

echo "Completed.... integron finder ------------------------------------------"

###############################################################################





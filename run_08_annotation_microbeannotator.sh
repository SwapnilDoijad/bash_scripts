#!/bin/bash
###############################################################################
# 08 annotation
###############################################################################
echo "started... step-08 annotation -------------------------------------------"
###############################################################################
## preliminary file preparation and directory creation
    source /home/swapnil/miniconda3/etc/profile.d/conda.sh
    conda activate microbeannotator

    echo "provide list file (for e.g. all)"
    echo "-------------------------------------------------------------------------"
    ls list.*.txt | sed 's/ /\n/g' | sed 's/list\.//g' | sed 's/\.txt//g'
    echo "-------------------------------------------------------------------------"
    read l
    list=$(echo "list.$l.txt")

    (mkdir results/08_annotation_microbeannotator) > /dev/null 2>&1
###############################################################################
for F1 in $(cat $list); do
echo results/08_annotation/raw_files/$F1/$F1.faa >> list.$l.MicrobeAnnotator.txt
done

/home/swapnil/tools/MicrobeAnnotator/MicrobeAnnotator \
-l list.$l.MicrobeAnnotator.txt \
-d /media/swapnil/databases/bioinfoDBs/MicrobeAnnotator_DB \
-o results/08_annotation_microbeannotator \
-m sword -p 4 -t 8 \
--kofam_bin /media/swapnil/databases/bioinfoDBs/kofamscan/kofam_scan-1.3.0 --light

conda deactivate
###############################################################################

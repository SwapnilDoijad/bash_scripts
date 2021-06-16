#!/bin/bash
###############################################################################
# 08 annotation
###############################################################################

echo "started... step-08 annotation -------------------------------------------"

###############################################################################
    source /home/swapnil/miniconda3/etc/profile.d/conda.sh
    conda activate myenv

    echo "provide list file (for e.g. all)"
    echo "---------------------------------------------------------------------"
    ls list.*.txt | awk -F'.' '{print $2}'
    echo "---------------------------------------------------------------------"
    read l
    list=$(echo "list.$l.txt")

    (mkdir results/08_annotation) > /dev/null 2>&1
    (mkdir results/08_annotation/raw_files) > /dev/null 2>&1

for F1 in $(cat $list); do

    echo "running.... step-08 annotation for..." $F1
    if [ -f results/08_annotation/raw_files/$F1/$F1.gbk ]; then
        echo "Annotation of $F1 alredy finished"
        else
        (mkdir results/08_annotation/raw_files/$F1) > /dev/null 2>&1
        (mkdir results/08_annotation/raw_files/$F1/tmp) > /dev/null 2>&1
                
        prokka --quiet --outdir results/08_annotation/raw_files/$F1 --force --prefix $F1 --addgenes --locustag $F1 --strain $F1 --rnammer results/04_assembly/all_fasta/"$F1".fasta

        V1=$(grep "bases:" results/08_annotation/raw_files/$F1/$F1.txt)
        V2=$(grep "rRNA:" results/08_annotation/raw_files/$F1/$F1.txt)
        V3=$(grep "tRNA:" results/08_annotation/raw_files/$F1/$F1.txt)
        V4=$(grep "gene:" results/08_annotation/raw_files/$F1/$F1.txt)
        V5=$(grep "CDS:" results/08_annotation/raw_files/$F1/$F1.txt)

        V6a=$(perl /home/swapnil/pipeline/tools/get_gc_content.pl results/08_annotation/raw_files/$F1/$F1.fna)
        V6=$(echo "GC:" $V6a)
        printf "$V1\n$V2\n$V3\n$V4\n$V5\n$V6" > results/08_annotation/raw_files/$F1/tmp/$F1.08_annotation.statistics.txt

        #------------------------------------------------------------------------------
        #gbk2csv
        (ugene --task=/home/swapnil/pipeline/tools/gbk2csv-Ugene.uwl --in=results/08_annotation/raw_files/$F1/$F1.gbk --out=results/08_annotation/raw_files/$F1/$F1.gbk.csv --format=csv) > /dev/null 2>&1
        awk -F'"' -v OFS='' '{ for (i=2; i<=NF; i+=2) gsub(",", "", $i) } 1' results/08_annotation/raw_files/$F1/$F1.gbk.csv | sed 's/"//g' > results/08_annotation/raw_files/$F1/$F1.gbk.internal_comma_quote_removed.csv
        #------------------------------------------------------------------------------

    fi
    echo "finished... step-08 annotation for..." $F1

done

conda deactivate
###############################################################################

echo "completed... step-08 annotation -----------------------------------------"

###############################################################################

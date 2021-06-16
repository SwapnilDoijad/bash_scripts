#!/bin/bash
###############################################################################
echo "started.... ANI-enveomics ----------------------------------------------"
###############################################################################
## initial inout and file preparation
    echo "provide list file (for e.g. all)"
    echo "---------------------------------------------------------------------"
    ls list.*.txt | sed 's/ /\n/g'
    echo "---------------------------------------------------------------------"
    read l
    list=$(echo "list.$l.txt")
    sfx=$(echo "_$l")

    echo "Genus to includes"
    echo "-------------------------------------------------------------------------"
    ls /media/swapnil/share/databases/TS_database/*.fasta | awk -F'/' '{print $NF}' | awk -F'_' '{print $1}' | sort -u
    echo "-------------------------------------------------------------------------"
    read genus
    TS_list=$(ls /media/swapnil/share/databases/TS_database/"$genus"*.fasta | awk -F'/' '{print $NF}' | sort -u | sed 's/\.fasta//g' )
    TS_path=$(echo "/media/swapnil/share/databases/TS_database")
    echo $TS_list > list.TS.txt

    (mkdir results/06_TW-ANI_enveomics_identification"$sfx") > /dev/null 2>&1
    (mkdir results/06_TW-ANI_enveomics_identification"$sfx"/raw_files) > /dev/null 2>&1
    (mkdir results/06_TW-ANI_enveomics_identification"$sfx"/raw_files/tmp) > /dev/null 2>&1

###############################################################################
    for TS in $(cat list.TS.txt); do
        (mkdir results/06_TW-ANI_enveomics_identification"$sfx"/raw_files/$TS) > /dev/null 2>&1
        if [ ! -f results/06_TW-ANI_enveomics_identification"$sfx"/raw_files/$TS.ANI.csv ]; then
            echo "subj quer two-way-ANI(%) TWA-fragm subj-fragm quer-fragm %Subj_covered %Quer_covered" > results/06_TW-ANI_enveomics_identification"$sfx"/raw_files/$TS.ANI.csv
        fi
        for F2 in $(cat $list); do
            if [ ! -f results/06_TW-ANI_enveomics_identification"$sfx"/raw_files/$TS/$TS.$F2.res.out_file ]; then
                echo "running $TS vs $F2"
                ( ruby /home/swapnil/tools/enveomics-master/Scripts/ani.rb --seq1 $TS_path/$TS.fasta --seq2 results/04_assembly/all_fasta/$F2.fasta --res results/06_TW-ANI_enveomics_identification"$sfx"/raw_files/$TS/$TS.$F2.res.out_file --quiet --threads 8 -p blast -w 1020 -s 150 -i 30 -l 714 ) > /dev/null 2>&1

                V1=$(awk 'FNR==3{print $3}' results/06_TW-ANI_enveomics_identification"$sfx"/raw_files/$TS/$TS.$F2.res.out_file | sed 's/\%//g')
                V2=$(awk 'FNR==3{print $7}' results/06_TW-ANI_enveomics_identification"$sfx"/raw_files/$TS/$TS.$F2.res.out_file)
                V3=$(awk 'FNR==1{print $8}' results/06_TW-ANI_enveomics_identification"$sfx"/raw_files/$TS/$TS.$F2.res.out_file)
                V4=$(awk 'FNR==2{print $8}' results/06_TW-ANI_enveomics_identification"$sfx"/raw_files/$TS/$TS.$F2.res.out_file)
                V5=$(echo 'scale=1;100*'$V2'/'$V3'' | bc )
                V6=$(echo 'scale=1;100*'$V2'/'$V4'' | bc )           
                echo "$TS $F2 $V1 $V2 $V3 $V4 $V5 $V6" >> results/06_TW-ANI_enveomics_identification"$sfx"/raw_files/$TS.ANI.csv
            else
            echo "ANI already calcualted for $TS and $F2"
            fi
        done
    done

(rm results/06_TW-ANI_enveomics_identification"$sfx"/tmp/results.csv) > /dev/null 2>&1
awk '$3>95 {print $0}' results/06_TW-ANI_enveomics_identification"$sfx"/raw_files/*.csv | sed '/\%Quer_covered/d' | awk '{print $1, $2, $3}' > results/06_TW-ANI_enveomics_identification"$sfx"/raw_files/tmp/results.csv

cat $list | sort > list1.tmp
cat results/06_TW-ANI_enveomics_identification"$sfx"/raw_files/tmp/results.csv | awk '{print $2}' | sort > list2.tmp
comm -3 list1.tmp list2.tmp | sed -e 's/\t//g' > results/06_TW-ANI_enveomics_identification"$sfx"/raw_files/tmp/unidentified.tmp
sed -i -e 's/^/unidentified\ /g' results/06_TW-ANI_enveomics_identification"$sfx"/raw_files/tmp/unidentified.tmp
cat results/06_TW-ANI_enveomics_identification"$sfx"/raw_files/tmp/results.csv results/06_TW-ANI_enveomics_identification"$sfx"/raw_files/tmp/unidentified.tmp > results/06_TW-ANI_enveomics_identification"$sfx"/results.csv
rm list1.tmp
rm list2.tmp
###############################################################################
echo "completed.... ANI-enveomics --------------------------------------------"
###############################################################################

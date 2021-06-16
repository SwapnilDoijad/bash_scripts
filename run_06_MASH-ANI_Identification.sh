#!/bin/bash
read -p "Hello Swapnil, There is a bug in the script, would you like to continue? Identified during Enteorobacter ncbi project. It appers that if isolate has >95% ANI to the two TS then, it will print either the first or seconf type strains. This might leas to miss-identification. The scrrep needs to modify to choose the best (highest ANI) TS"
###############################################################################
#06 ANI_Ident
## to-do: Create a central sketch-datbase for type-strains (remember sketch takes time, distance calculation do not)
###############################################################################
echo "started.... step-06 ANI_Ident --------------------------------------------"
###############################################################################
echo "provide list file (for e.g. all)"
echo "-------------------------------------------------------------------------"
ls list.*.txt | sed 's/ /\n/g'
echo "-------------------------------------------------------------------------"
read l
list=$(echo "list.$l.txt")

echo "Genus to includes"
echo "-------------------------------------------------------------------------"
ls /media/swapnil/share/databases/TS_database/*.fasta | awk -F'/' '{print $NF}' | awk -F'_' '{print $1}' | sort -u
echo "-------------------------------------------------------------------------"
read genus
TS_list=$(ls /media/swapnil/share/databases/TS_database/"$genus"*.fasta | awk -F'/' '{print $NF}' | sort -u | sed 's/\.fasta//g' )
TS_path=$(echo "/media/swapnil/share/databases/TS_database")

(mkdir results/06_ANI_Ident) > /dev/null 2>&1
(mkdir results/06_ANI_Ident/tmp) > /dev/null 2>&1
(mkdir results/06_ANI_Ident/fasta) > /dev/null 2>&1
(mkdir results/06_ANI_Ident/msh) > /dev/null 2>&1
(mkdir results/06_ANI_Ident/msh/distance) > /dev/null 2>&1
(mkdir results/06_ANI_Ident/msh/distance/tmp) > /dev/null 2>&1
###############################################################################
## genome-sketch type strains (if not prensted earlier)
    for ref_TS in $TS_list; do
    if [ ! -f /media/swapnil/share/databases/TS_database/msh/"$ref_TS".msh ]; then
    echo sketching "$ref_TS" 
    (mash sketch -o $TS_path/msh/"$ref_TS".msh $TS_path/"$ref_TS".fasta)> /dev/null 2>&1
    fi
    done
## genome-sketch for test isolates
    for F1 in $(cat $list); do
    cp results/04_assembly/all_fasta/"$F1".fasta results/06_ANI_Ident/fasta/
    echo sketching "$F1" 
    (mash sketch -o results/06_ANI_Ident/msh/"$F1".msh results/06_ANI_Ident/fasta/"$F1".fasta)> /dev/null 2>&1
    done
## mash distance for type strains 
    for ref_TS in $TS_list; do
    echo running mash-distance "$ref_TS"
    mash dist $TS_path/msh/"$ref_TS".msh results/06_ANI_Ident/msh/*.msh > results/06_ANI_Ident/msh/distance/"$ref_TS".mash-distance.tab
    sed -i "s/results\/06_ANI_Ident\/fasta\///g" results/06_ANI_Ident/msh/distance/"$ref_TS".mash-distance.tab
    sed -i "s/\.fasta//g" results/06_ANI_Ident/msh/distance/"$ref_TS".mash-distance.tab
    done
## mash distance for test isolates
    (rm results/06_ANI_Ident/tmp/results.csv)> /dev/null 2>&1
    for ref_TS in $TS_list; do
    for F1 in $(cat $list) ; do
    ANI=$(awk '$2=="'$F1'" {print (1-$3)*100}' results/06_ANI_Ident/msh/distance/"$ref_TS".mash-distance.tab | awk -F'.' '{print $1}')
    if [ "$ANI" -ge "95" ] ; then
    ANI2=$(awk '$2=="'$F1'" {print (1-$3)*100}' results/06_ANI_Ident/msh/distance/"$ref_TS".mash-distance.tab)
    echo "$F1 $ref_TS $ANI2" >> results/06_ANI_Ident/tmp/results.csv
    fi
    done
    done
## summarise results
    echo "summarising..."
    awk '{print $1}' results/06_ANI_Ident/tmp/results.csv > results/06_ANI_Ident/tmp/tmp.1.tmp
    cat $list > results/06_ANI_Ident/tmp/tmp.2.tmp
    cat results/06_ANI_Ident/tmp/tmp.1.tmp results/06_ANI_Ident/tmp/tmp.2.tmp > results/06_ANI_Ident/tmp/tmp.3.tmp
    awk '{!seen[$0]++};END{for(i in seen) if(seen[i]==1)print i}' results/06_ANI_Ident/tmp/tmp.3.tmp | sed 's/$/ unidentified ????/' >> results/06_ANI_Ident/tmp/results.csv

    (rm results/06_ANI_Ident/tmp/results.2.csv)> /dev/null 2>&1
    for F1 in $(cat $list); do
    grep "$F1" results/06_ANI_Ident/tmp/results.csv | sort -k3,3r | awk 'FNR==1 {print $1, $2, $3}' >> results/06_ANI_Ident/tmp/results.2.csv
    done

    (rm results/06_ANI_Ident/tmp/results.stat.csv)> /dev/null 2>&1
    total_isolates=$(cat $list | wc -l )
    echo $TS_list > results/06_ANI_Ident/tmp/list.tmp
    echo "unidentified" >> results/06_ANI_Ident/tmp/list.tmp
    for ref_TS in $(cat results/06_ANI_Ident/tmp/list.tmp); do
    V1=$(grep "$ref_TS" results/06_ANI_Ident/tmp/results.2.csv | wc -l)
    V2=$(bc -l <<< 'scale=2; '$V1'*100/'$total_isolates'')
    echo $ref_TS $V1 $V2 >> results/06_ANI_Ident/tmp/results.stat.csv
    done
    rm results/06_ANI_Ident/tmp/list.tmp

    sort -k2,2rn results/06_ANI_Ident/tmp/results.stat.csv > results/06_ANI_Ident/tmp/results.stat.2.csv
    sed -i '1 i\Isolate Species ANI%' results/06_ANI_Ident/tmp/results.2.csv
    sed -i '1 i\Species Total Percentage' results/06_ANI_Ident/tmp/results.stat.csv

    ssconvert results/06_ANI_Ident/tmp/results.2.csv results/06_ANI_Ident/results.xlsx
    ssconvert results/06_ANI_Ident/tmp/results.stat.2.csv results/06_ANI_Ident/results.stat.xlsx

    awk -F'_' '{print $1,$2}' results/06_ANI_Ident/tmp/results.csv | awk '{print $1,$1,$2,$3}' | sed "s/ /_/3" | sed "s/ /(/2" | sed -e "s/$/)/" > results/06_ANI_Ident/renaming.txt
    sed -i '1 i\#old\tnew' results/06_ANI_Ident/renaming.txt
    sed -i 's/ /\t/g' results/06_ANI_Ident/renaming.txt
###############################################################################
echo "completed.... step-06 ANI_Ident ------------------------------------------"
###############################################################################
